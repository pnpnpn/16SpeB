#!/usr/bin/perl -w
use strict;
use diagnostics;

#Global constants
use constant DEBUG0 => 1;
#use lib "$ENV{PROJECT_DIRECTORY}/../uri/bacteria/r16s/r16s_scripts/";

#libraries
use Benchmark qw(:all);
use File::Basename;
use File::Temp;
use FastaIO;
use FastaFilter;
use MyRandom;
use MyMath;
use R16sHelper;
use List::Util qw(sum min max);
use POSIX qw(ceil floor);

my $TMP_DIR = '/tmp/';
my $ALIGN_EXE = "palign/palign.out";

if(!-f $ALIGN_EXE) {
	die "Missing palign.out file. Try running ./compile_palign";
}

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
usage: <program> --dir=<DIR> [OPTIONS]

Compute the intra-species PID distribution.

OPTIONS:
--iters=<INT>        Number of random iterations for each species when sampling
--exact=<INT>        Threshold of using exact method before random sampling 
--randseed=<INT>     Random seed
--genus=<STRING>     Genus-of-interest

PID extraction (choose one):
--aggregate          aggregate all the PIDs 

./compute_intra_species_pid_distrib.pl --iters=100 --dir=test_Mycoplasma_hominis --genus=Mycoplasma

Created on 6/16/2010

Example: 

USAGE
		exit(1);
	}

	my $user_randseed = undef;
	my $fsa_fname = undef;
	my $genus = undef;
	my $numseqs_before_rnd = 999999;
	my $rnd_iters = undef;
	my $fsadir = undef;
	my $is_aggregate = undef;
	foreach my $a (@ARGV) {
		if( $a =~ /^--randseed=(\d+)/) {
			$user_randseed = $1;
		}
		elsif( $a =~ /^--dir=(\S+)/) {
			$fsadir=$1;
		}
		elsif( $a =~ /^--genus=(\S+)/) {
			$genus=$1;
		}
		elsif( $a =~ /^--exact=(\d+)/) {
			$numseqs_before_rnd = $1;
		}
		elsif( $a =~ /^--iters=(\d+)/) {
			$rnd_iters = $1;
		}
		elsif( $a =~ /^--aggregate/) {
			$is_aggregate = 1;
		}
		else {
			die "Unrecognized parameter: $a\n";
		}
	}

	if(!defined($genus)) {
		die "genus must be defined";
	}
	if(!defined($fsadir)) {
		die "fsadir must be defined";
	}

	if(!-d $fsadir) {
		die "Error: $fsadir does not exist";
	}

	#if a new random seed is read
	if(defined($user_randseed)) {
		$RAND_SEED = $user_randseed;
	}
	$RAND_MT_GEN->srand($RAND_SEED);
	print STDERR "Random seed: " . $RAND_SEED . "\n";

	opendir(my $dh, $fsadir) || die "Cannot open directory $fsadir: $!";
	my @fsalst = ();
	foreach my $f (readdir($dh)) {
		if($f =~ /\.fsa$/ ) {
			push(@fsalst, $f);
		}
	}
	closedir($dh);

	my $count = 0;
	my %pids_hash= ();
	my %numseqs_hash = ();
	my @aggregate_pids = ();
	foreach my $f (@fsalst) {
		if($f =~ /^${genus}_(\w+)\.fsa/) {
			my $species = $1;
			my $mode = undef;
			my $fsa = parseFa("$fsadir/$f");
			my $numseqs = scalar(@$fsa);
			if($numseqs > $numseqs_before_rnd) {
				$mode = "-rand-pair $rnd_iters";
				print STDERR "WARNING: using --rand-pair mode\n";
			}
			else {
				$mode = "-all-pair";
			}
			my $outfname = sprintf("%s_%s.txt", $genus, $species);
			my $cmd = "$ALIGN_EXE $fsadir/$f $mode -s $RAND_SEED -quiet 1>$fsadir/$outfname 2>&1";
			print STDERR "$cmd\n\n";
			system("$cmd");
			
			#pid list
			my $pid_over_alignlen_fn = sprintf("%s_%s_pid_over_alignlen.txt", $genus, $species);
			$cmd = "grep 'PID over align' $fsadir/$outfname | awk {'print \$4'} > $fsadir/$pid_over_alignlen_fn";
			print STDERR "$cmd\n\n";
			system($cmd);

			#parse pid file
			my $sorted = extract_and_sort_pids("$fsadir/$pid_over_alignlen_fn");

			push(@aggregate_pids, @$sorted);
			$numseqs_hash{$species} = $numseqs;

			if(!$is_aggregate) {
				#minimum
				$pids_hash{$species} = $sorted;
			}

			#count
			$count++;
		}
		else {
			die "Invalid FASTA file found: $f\n";
		}
	}


	print STDERR "\n";
	print STDERR "Number of unique species: $count\n";
	print STDERR "\n";



	my %minpid_hash = ();
	print STDERR "PIDs (min, 5%-quantile, median, 95%-quantile, max):\n";
	foreach my $species (sort keys %pids_hash){
		my $sorted = $pids_hash{$species};
		my $minpid = min(@$sorted);
		my $quantile05 = compute_empirical_quantile($sorted, 0.05);
		my $quantile95 = compute_empirical_quantile($sorted, 0.95);
		my $median = compute_empirical_quantile($sorted, 0.50);
		my $maxpid = max(@$sorted);
		$minpid_hash{$species} = $minpid;
		printf STDERR (
			"(%.5lf, %.5lf, %.5lf, %.5lf, %.5lf)  $genus $species with %d sequences (%d pairs)\n",
			$minpid,
			$quantile05,
			$median, 
			$quantile95,
			$maxpid,
            $numseqs_hash{$species},
			scalar(@$sorted),
		);

	}
	print STDERR "\n";

	if($is_aggregate) {
		printf join("\n", @aggregate_pids) . "\n";
	}
	else {
		foreach my $species (sort keys %minpid_hash){
			printf $minpid_hash{$species} . "\n";
		}
	}
	
}

sub extract_and_sort_pids
{
	my ($fname) = @_;

	my @lst = ();
	open(IN, "<$fname") or die "Cannot open $fname for read: $!";
	foreach my $ln (<IN>) {
		$ln =~ s/\s//g;
		if($ln =~ /^[\d\.]+$/) {
			push(@lst, $ln);
		}
		elsif(length($ln) == 0) {
			#skip
		}
		else {
			die "Cannot parse line: $ln";
		}
	}
	close(IN);
	my @sorted = sort {$a <=> $b} @lst;
	return \@sorted;
}



