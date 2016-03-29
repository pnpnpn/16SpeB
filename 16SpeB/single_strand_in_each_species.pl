#!/usr/bin/perl -w
use strict;
use diagnostics;

#Global constants
use constant DEBUG0 => 1;

#libraries
use Benchmark qw(:all);
use File::Basename;
use File::Temp;
use FastaIO;
use CnsgSampling;
use MyRandom;
use POSIX qw(ceil floor);

my $TMP_DIR = '/tmp/';

my $GAPCHAR = 4;
my $NUMALPHAS = 4;

my $RAND_SEED = time ^ $$ ^ unpack "%L*", `ps axww | gzip`;
my $RAND_MT_GEN = MyRandom->new();

my $start_benchmark = new Benchmark();
main();
my $end_benchmark = new Benchmark();

print STDERR "\n";
printf STDERR ("Benchmark: %s\n", timestr(timediff($end_benchmark, $start_benchmark), 'all'));


sub main 
{
	print STDERR `date` . "\n";
	if(scalar(@ARGV) == 0) {
		print "usage: <program> --fsa=<FSA>\n";
		print "\n";
		print "Randomly keep a single strand in each species (specifically for 16s downloaded from rdp).\n";
		print "\n";
		print "OPTIONS:\n";
		print "--randseed=<INT>\n";
		print "--genus=<STRING>\n";
		print "\n";
		print "Created on 5/19/2010\n";
		print "\n";
		exit(1);
	}


	my $user_randseed = undef;
	my $fsa_fname = undef;
	my $genus = undef;
	foreach my $a (@ARGV) {
		if( $a =~ /^--randseed=(\d+)/) {
			$user_randseed = $1;
		}
		elsif( $a =~ /^--fsa=(\S+)/) {
			$fsa_fname=$1;
		}
		elsif( $a =~ /^--genus=(\S+)/) {
			$genus=$1;
		}
		else {
			die "Unrecognized parameter: $a\n";
		}
	}
	
	if(!defined($genus)) {
		die "genus must be defined";
	}
	#if a new random seed is read
	if(defined($user_randseed)) {
		$RAND_SEED = $user_randseed;
	}
	$RAND_MT_GEN->srand($RAND_SEED);
	print STDERR "Random seed: " . $RAND_SEED . "\n";

	my $fsa = parseFa($fsa_fname);
	my $numseqs = scalar(@$fsa);

	print STDERR "CMD: $0 " . join(' ', @ARGV) . "\n";
	print STDERR "Number of sequences: $numseqs\n";
	print STDERR "\n";

	$fsa = filter_by_genus(*STDERR, $fsa, $genus);

	my $new_fsa = keep_single_strand_in_each_species(*STDERR, $fsa);
	writeFa(*STDOUT, $new_fsa);

	print STDERR "\n";
	print STDERR "Number of unique species: " . scalar(@$new_fsa)."\n";
}

sub filter_by_genus
{
	my ($fptr, $fsa, $genus) = @_;

	my @new_fsa = ();
	foreach my $s (@$fsa) {
		if($s->{tag} =~ /$genus/i) {
			push(@new_fsa, $s);
		}
		else {
			print $fptr ("Removed unidentified: ". $s->{tag} . "\n");
		}
	}
	print $fptr "\n";
	return \@new_fsa;
	
}

sub keep_single_strand_in_each_species 
{
	#fptr is only for outputing "stderr" output. The actual FASTA
	#file is returned by this function.

	my ($fptr, $fsa) = @_;
	my %species_bins = ();
	$species_bins{'uncultured'} = [];

	my $genus0 = undef;
	foreach my $s (@$fsa) {
		if($s->{tag} =~ /uncultured/) {
			push(@{$species_bins{'uncultured'}}, $s);
		}
		elsif(
			$s->{tag} =~ /^S\d+ (\w+) (\w+)[\s\;]/
			|| $s->{tag} =~ /^S\d+ (\w+) sp\. ([\d\w\-]+)[\s\;]/
		){
			my $genus = $1;
			my $species = $2;
			if(!exists($species_bins{$species})) {
				$species_bins{$species} = [];
			}
			push(@{$species_bins{$species}}, $s);

			#sanity check
			if(!defined($genus0)) {
				$genus0 = $genus;
			}
			if($genus0 ne $genus) {
				die "Unequal genus found $genus0 != $genus";
			}
		}
		else {
			print STDERR "WARNING: Cannot parse tag: " . $s->{tag} . "\n";
		}
	}
	print STDERR "\n";

	my @species_names = sort keys %species_bins;
	my @new_fsa = ();
	foreach my $name (@species_names) {
		my $num_strands = scalar(@{$species_bins{$name}});
		my @taglst = ();
		foreach my $s (@{$species_bins{$name}}) {
			my $tag = $s->{tag};
			$tag =~ s/\s/\_/g;
			push(@taglst, $tag);
		}
		printf $fptr ("%3d strands found for species $name: %s\n",
			$num_strands,
			join(" ", @taglst),
		);

		if($num_strands > 0) {
			my $rnd_ind = floor($RAND_MT_GEN->rand() * $num_strands);
			push(@new_fsa, $species_bins{$name}->[$rnd_ind]);
		}
	}

	return \@new_fsa;
}


