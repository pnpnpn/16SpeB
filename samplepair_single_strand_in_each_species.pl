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
use FastaFilter;
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
		print <<USAGE;
usage: <program> --fsa=<FSA>

Randomly sample a pair of strands from two different species (specifically 
for 16s downloaded from rdp). It outputs a FASTA file that lists the 
pairs.

OPTIONS:
--numpairs=<INT>     Number of pairs to sample
--randseed=<INT>     Random seed
--genus=<STRING>     Genus-of-interest

Created on 6/02/2010

USAGE
		exit(1);
	}

	my $user_randseed = undef;
	my $fsa_fname = undef;
	my $genus = undef;
	my $numpairs= undef;
	foreach my $a (@ARGV) {
		if( $a =~ /^--randseed=(\d+)/) {
			$user_randseed = $1;
		}
		elsif( $a =~ /^--fsa=(\S+)/) {
			$fsa_fname=$1;
		}
		elsif( $a =~ /^--numpairs=(\d+)/) {
			$numpairs=$1;
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
	print STDERR "Number of pairs: $numpairs\n";
	print STDERR "\n";

	$fsa = filter_by_genus(*STDERR, $fsa, $genus);

	my @badwords = (
		' sp\.? ', 
		' genomosp\.? ', 
		' uncultured ', 
		' Persistence ', 
		' stable enrichment ', 
	);
	$fsa = filter_out_keywords(*STDERR, $fsa, \@badwords);

	my $species_hash = bin_strands_by_species(*STDERR, $fsa);

	my $new_fsa = sample_pairs_from_each_species(*STDERR, $species_hash, $numpairs);
	writeFa(*STDOUT, $new_fsa);

	print STDERR "\n";
	print STDERR "Number of unique species: " . scalar(keys %$species_hash)."\n";
	print STDERR "\n";
	
	my $sum = 0;
	my $minlen = length($new_fsa->[0]->{seq}) ;
	my $maxlen = 0;
	foreach my $s (@$new_fsa) {
		my $len =  length($s->{seq});
		$sum += $len;
		if($minlen > $len) {
			$minlen = $len;
		}
		if($maxlen < $len) {
			$maxlen = $len;
		}
	}
	printf STDERR ("Min sequence length: %d\n", $minlen);
	printf STDERR ("Average sequence length: %d\n", $sum / scalar(@$new_fsa));
	printf STDERR ("Max sequence length: %d\n", $maxlen);
}

sub sample_pairs_from_each_species 
{
	my ($fptr, $species_hash, $numpairs) = @_;

	my @new_fsa = ();
	my @species_names = sort keys %$species_hash;
	for(my $i = 0; $i < $numpairs; $i++) {
		my $rnd_ind1 = undef;
		my $rnd_ind2 = undef;
		do {
			$rnd_ind1 = floor($RAND_MT_GEN->rand() * scalar(@species_names));
			$rnd_ind2 = floor($RAND_MT_GEN->rand() * scalar(@species_names));
		} while($rnd_ind1 == $rnd_ind2);

		my $strand_lst1 = $species_hash->{$species_names[$rnd_ind1]};
		my $strand_lst2 = $species_hash->{$species_names[$rnd_ind2]};

		my $strand1 = sample_strand_from_list($strand_lst1);
		my $strand2 = sample_strand_from_list($strand_lst2);

		push(@new_fsa, $strand1);
		push(@new_fsa, $strand2);
	}
	return \@new_fsa;
}

sub sample_strand_from_list 
{
	my ($lst) = @_;
	if(scalar(@$lst) == 0) {
		die "Invalid number of strands found";
	}
	my $rnd_ind = floor($RAND_MT_GEN->rand() * scalar(@$lst));
	return $lst->[$rnd_ind];
}


sub bin_strands_by_species
{
	#fptr is only for outputing "stderr" output. 

	my ($fptr, $fsa) = @_;
	my %species_bins = ();

	my $genus0 = undef;
	foreach my $s (@$fsa) {
		$s->{tag} =~ s/\"/ /g;
		my $tag = $s->{tag};
		$tag =~ s/ subsp\. / /g;
		$tag .= ' ';
		if(
			$tag =~ /^S\d+\s+(\w+) (\w+)[\s\;]/
			|| $tag =~ /^silva\|\S+\s+(\w+) (\w+)[\s\;]/
			|| $tag =~ /^\d+ \S+\s+(\w+) (\w+)[\s\;]/
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
			die;
		}
	}
	print STDERR "\n";

	my @species_names = sort keys %species_bins;
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
	}
	return \%species_bins;
}




