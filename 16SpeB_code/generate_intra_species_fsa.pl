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
	print STDERR "$0 " . join(' ', @ARGV) . "\n\n";
	if(scalar(@ARGV) == 0) {
		print <<USAGE;
usage: <program> --fsa=<FSA> --outdir=<DIR>

Generate an individual FASTA file for each species. This is useful for
analyzing intra-species variation.

OPTIONS:
--min-strands                Minimum number of strands allowed for adding to FASTA
--randseed=<INT>             Random seed
--genus=<STRING>             Genus-of-interest
--bi-nomenclature=<STRING>   species-of-interest (e.g. Saccharomyces_cerevisiae)

Created on 6/14/2010

USAGE
		exit(1);
	}

	my $user_randseed = undef;
	my $fsa_fname = undef;
	my $genus = undef;
	my $min_strands_threshold = 0;
	my $outdir = undef;
	my $bi_nomenclature = undef;
	foreach my $a (@ARGV) {
		if( $a =~ /^--randseed=(\d+)/) {
			$user_randseed = $1;
		}
		elsif( $a =~ /^--fsa=(\S+)/) {
			$fsa_fname=$1;
		}
		elsif( $a =~ /^--outdir=(\S+)/) {
			$outdir=$1;
		}
		elsif( $a =~ /^--genus=(\S+)/) {
			$genus=$1;
		}
		elsif( $a =~ /^--min-strands=(\d+)/) {
			$min_strands_threshold=$1;
		}
		elsif( $a =~ /^--bi-nomenclature=(\S+)/) {
			$bi_nomenclature=$1;
		}
		else {
			die "Unrecognized parameter: $a\n";
		}
	}

	if(!(defined($genus) xor defined($bi_nomenclature))) {
		die "--genus and --bi-nomenclature can not defined together";
	}

	my $species_interest = undef;
	($genus, $species_interest) = split('_', $bi_nomenclature);

	if(!defined($genus)) {
		die "genus must be defined";
	}

	if(!defined($outdir)) {
		die "outdir must be defined";
	}

	if(!-d $outdir) {
		print STDERR ("Creating directory $outdir...\n");
		mkdir($outdir) or die "cannot create directory";
	}
	if(!-d $outdir) {
		die "Error: $outdir does not exist";
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
	print STDERR "Minimum number of strands threshold: $min_strands_threshold\n";
	print STDERR "\n";

	my $species_hash = undef;
	if(defined($bi_nomenclature)) {
		$species_hash = {};
		$species_hash->{$species_interest} = $fsa;
	}
	else {
		$species_hash = bin_strands_by_species($fsa);
	}
	print_species_bins(*STDERR, $species_hash);
	print STDERR "\n";

	foreach my $species (keys %$species_hash) {
		my $numstrands = scalar(@{$species_hash->{$species}});
		if($numstrands >= $min_strands_threshold) {
			my $fname = sprintf("%s_%s.fsa",$genus,$species);
			open(my $outfh, "> $outdir/$fname") or die "cannot open $outdir/$fname for write: $!";
			writeFa($outfh, $species_hash->{$species});
			close($outfh);
		}
		else {
			print STDERR ("Species $species removed because it only have $numstrands strands\n");
		}
	}

	print STDERR "\n";
	print STDERR "Number of unique species: " . scalar(keys %$species_hash)."\n";
	print STDERR "\n";

	my @lenlst = ();
	foreach my $species (keys %$species_hash) {
		foreach my $s (@{$species_hash->{$species}}) {
			push(@lenlst, length($s->{seq}));
		}
	}
	@lenlst = sort {$a<=> $b}@lenlst;
	printf STDERR ("Min sequence length: %d\n", $lenlst[0]);
	printf STDERR ("5%%-quantile sequence length: %d\n", compute_empirical_quantile(\@lenlst, 0.05));
	printf STDERR ("Median sequence length: %d\n", compute_empirical_quantile(\@lenlst, 0.5));
	printf STDERR ("95%%-quantile sequence length: %d\n", compute_empirical_quantile(\@lenlst, 0.95));
	printf STDERR ("Max sequence length: %d\n", $lenlst[scalar(@lenlst)-1]);
}
