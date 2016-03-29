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
use MyMath;
use R16sHelper;
use List::Util qw(sum);
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
	print STDERR "$0 " . join(" ", @ARGV) . "\n\n";
	if(scalar(@ARGV) == 0) {
		print <<USAGE;
usage: <program> --fsa=<FSA>

Filter duplicates within intra-species

OPTIONS:
--randseed=<INT>     random seed
--genus=<STRING>     genus-of-interest (note this will filter by genus)
--species=<STRING>   species-of-interest (doesn't do any type of filtering and only looks at one bin)
--auxin=<FILE>       auxiliary file to filter (numseqs/seqlen in same order)
--auxout=<FILE>      output after filtering auxiliary file
--nobadwords         does not filter 'badwords'

Example:
./filter_duplicates_within_intra_species.pl --fsa=test_Mycoplasma_hominis/Mycoplasma_hominis_raw_16s.fsa --species=Mycoplasma_hominis --nobadwords --auxin=test_Mycoplasma_hominis/Mycoplasma_hominis_raw_v6.fsa --auxout=test_Mycoplasma_hominis/Mycoplasma_hominis_nondup_v6.fsa

Created on 6/14/2010

USAGE
		exit(1);
	}

	my $user_randseed = undef;
	my $fsa_fname = undef;
	my $genus = undef;
	my $species_interest= undef;
	my $use_badwordsfilter = 1;
	my $auxout_fname= undef;
	my $auxin_fname= undef;
	foreach my $a (@ARGV) {
		if( $a =~ /^--randseed=(\d+)/) {
			$user_randseed = $1;
		}
		elsif( $a =~ /^--fsa=(\S+)/) {
			$fsa_fname=$1;
		}
		elsif( $a =~ /^--auxin=(\S+)/) {
			$auxin_fname=$1;
		}
		elsif( $a =~ /^--auxout=(\S+)/) {
			$auxout_fname=$1;
		}
		elsif( $a =~ /^--genus=(\S+)/) {
			$genus=$1;
		}
		elsif( $a =~ /^--nobadwords/) {
			$use_badwordsfilter = 0;
		}
		elsif( $a =~ /^--species=(\S+)/) {
			$species_interest=$1;
		}
		else {
			die "Unrecognized parameter: $a\n";
		}
	}

	if(!(defined($species_interest) xor defined($genus))) {
		die "either --genus or --species should be chosen";
	}

	if(!defined($genus)) {
		#die "genus must be defined";
		#print STDERR "WARNING: genus is not defined\n";
	}

	if(defined($auxin_fname) xor defined($auxout_fname)) {
		die "Must specify both auxin and auxout.";
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

	if(defined($genus)) {
		$fsa = filter_by_genus(*STDERR, $fsa, $genus);
	}

	my @badwords = (
		' sp[\s\;\.]',
		' genomosp[\s\;\.]',
		' uncultured ',
		' Persistence ',
		' stable enrichment ',
	);
	if($use_badwordsfilter) {
		$fsa = filter_out_keywords(*STDERR, $fsa, \@badwords);
	}

	my $species_hash = undef;
	if(defined($species_interest)) {
		$species_hash = {};
		$species_hash->{$species_interest} = $fsa;
	}
	elsif(defined($genus)) {
   		$species_hash = bin_strands_by_species($fsa);
	}
	else {
		die "invalid species_hash not set";
	}

	print STDERR "Before filtering:\n";
	print_species_bins(*STDERR, $species_hash);
	print STDERR "\n";
	$species_hash = filter_duplicates_within_species_bins($species_hash);
	print STDERR "\n";
	print STDERR "After filtering:\n";
	print_species_bins(*STDERR, $species_hash);
	print STDERR "\n";

	print STDERR "\n";
	print STDERR "Number of unique species: " . scalar(keys %$species_hash)."\n";
	print STDERR "\n";

	my @newfsa = ();
	my @lenlst = ();
	foreach my $species (sort keys %$species_hash) {
		foreach my $s (@{$species_hash->{$species}}) {
			push(@newfsa, $s);
			push(@lenlst, length($s->{seq}));
		}
	}
	writeFa(*STDOUT, \@newfsa);

	#auxiliary file
	if(defined($auxin_fname)) {
		my $auxfsa= parseFa($auxin_fname);
		if(scalar(@$auxfsa) != $numseqs) {
			print
			die "Unequal number of sequences for auxin and input";
		}
		$auxfsa = filter_auxin_for_duplicates(\@newfsa, $auxfsa);
		open(my $auxfh, "> $auxout_fname") or die "Cannot write to $auxout_fname: $!";
		writeFa($auxfh, $auxfsa);
		close($auxfh);
	}

	@lenlst = sort {$a <=> $b} @lenlst;
	printf STDERR ("Min sequence length: %d\n", $lenlst[0]);
	printf STDERR ("5%%-quantile sequence length: %d\n", compute_empirical_quantile(\@lenlst, 0.05));
	printf STDERR ("Median sequence length: %d\n", compute_empirical_quantile(\@lenlst, 0.5));
	printf STDERR ("95%%-quantile sequence length: %d\n", compute_empirical_quantile(\@lenlst, 0.95));
	printf STDERR ("Max sequence length: %d\n", $lenlst[scalar(@lenlst)-1]);
}

sub filter_auxin_for_duplicates
{
	my ($nondupfsa, $auxfsa) = @_;

	my %hash = ();
	foreach my $s (@$auxfsa) {
		if(exists($hash{$s->{tag}})) {
			die "Duplicate tags in auxiliary file: " . $s->{tag} . "\n";
		}
		$hash{$s->{tag}} = $s;
	}

	my @new_auxfsa = ();
	foreach my $s (@$nondupfsa) {
		if(exists($hash{$s->{tag}})) {
			push(@new_auxfsa, $hash{$s->{tag}});
		}
	}
	return \@new_auxfsa;
}
