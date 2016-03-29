#!/usr/bin/perl -w
use strict;
use diagnostics;

use lib "$ENV{PROJECT_DIRECTORY}/../uri/bacteria/r16s/r16s_scripts/";

#Global constants
use constant DEBUG0 => 1;

#libraries
use Benchmark qw(:all);
use File::Basename;
use File::Temp;
use FastaIO;
use FastaFilter;
use ConsensusSeq;
use MyRandom;
use MyMath;
use MatrixPermutation;
use PwmDistance;
use POSIX qw(ceil floor);
use List::Util qw(sum min max);
use R16sHelper;

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
usage: <program> --fsa=<FSA> --leftstr=<STRING> --midstr=<STRING> --rightstr=<STRING> 

Isolates a region where endpoints are best matches of 2 user-specified k-mers.

OPTIONS:
leftstr                  left k-mer
midstr                   mid k-mer
rightstr                 right k-mer

leftcoord                literature left coordinate
midcoord                 literature mid coordinate
rightcoord               literature right coordinate
mismatch                 maximum number of mismatch per k-mer (default: 2)
trimpos=<INT>            count len from left of mid and then truncate
trimneg=<INT>            count len from right of mid and then truncate

offsetleft=<INT>         number of offset positions allowed for leftstr
offsetright=<INT>        number of offset positions allowed for rightstr

longlen-fsa=<FILE>       writes to file the filtered seqs with of longer length
longlen-trimpos=<INT>    trim in the positive direction based on the mid
longlen-trimneg=<INT>    trim in the negative direction based on the mid

example: r16s_scripts/isolate_region_by_multimarkers.pl --fsa=Lactobacillus_seqpair.fsa 

Created on 6/22/2010

USAGE
		exit(1);
	}

	my $fsa_fname = undef;

	my $leftmark = create_nullmarker(); 
	my $midmark = create_nullmarker();
	my $rightmark = create_nullmarker();
	#str, coord, arr, isValid
	

	my $params = {
		max_mismatch => 3,
		trimpos => undef,
		trimneg => undef,
		offsetleft => undef,
		offsetright => undef,

		longlen_fname => undef,
		longlen_trimpos => undef,
		longlen_trimneg => undef,

		leftmark => $leftmark,
		midmark => $midmark, 
		rightmark => $rightmark,

	};
	foreach my $a (@ARGV) {
		if( $a =~ /^--leftstr=(\S+)/) {
			$params->{leftmark}->{str} = $1;
			$params->{leftmark}->{isValid} = 1;
		}
		elsif( $a =~ /^--rightstr=(\S+)/) {
			$params->{rightmark}->{str} = $1;
			$params->{rightmark}->{isValid} = 1;
		}
		elsif( $a =~ /^--midstr=(\S+)/) {
			$params->{midmark}->{str} = $1;
			$params->{midmark}->{isValid} = 1;
		}
		elsif( $a =~ /^--leftcoord=(\d+)/) {
			$params->{leftmark}->{coord} = $1;
			$params->{leftmark}->{isValid} = 1;
		}
		elsif( $a =~ /^--rightcoord=(\d+)/) {
			$params->{rightmark}->{coord} = $1;
			$params->{rightmark}->{isValid} = 1;
		}
		elsif( $a =~ /^--midcoord=(\d+)/) {
			$params->{midmark}->{coord} = $1;
			$params->{midmark}->{isValid} = 1;
		}
		elsif( $a =~ /^--fsa=(\S+)/) {
			$fsa_fname=$1;
		}
		elsif( $a =~ /^--mismatch=(\d+)/) {
			$params->{max_mismatch} = $1;
		}
		elsif( $a =~ /^--trimneg=(\d+)/) {
			$params->{trimneg} = $1;
		}
		elsif( $a =~ /^--trimpos=(\d+)/) {
			$params->{trimpos} = $1;
		}
		elsif( $a =~ /^--offsetleft=(\d+)/) {
			$params->{offsetleft} = $1;
		}
		elsif( $a =~ /^--offsetright=(\d+)/) {
			$params->{offsetright} = $1;
		}
		elsif( $a =~ /^--longlen-fsa=(\S+)/) {
			$params->{longlen_fname} = $1;
		}
		elsif( $a =~ /^--longlen-trimneg=(\d+)/) {
			$params->{longlen_trimneg} = $1;
		}
		elsif( $a =~ /^--longlen-trimpos=(\d+)/) {
			$params->{longlen_trimpos} = $1;
		}
		else {
			die "Unrecognized parameter: $a\n";
		}
	}
	if(defined($params->{longlen_trimneg}) xor defined($params->{longlen_trimpos})) {
		die "Both longlen-trim pos/neg must be defined";
	}
	if(defined($params->{trimneg}) xor defined($params->{trimpos})) {
		die "Both trim pos/neg must be defined";
	}
	if($params->{longlen_trimneg} < 0 || $params->{longlen_trimpos} < 0) {
		die "trim pos/neg must be positive";
	}
	
	#if a new random seed is read
	#if(defined($user_randseed)) {
	#	$RAND_SEED = $user_randseed;
	#}
	$RAND_MT_GEN->srand($RAND_SEED);
	#print STDERR "Random seed: " . $RAND_SEED . "\n";

	my $fsa = parseFa($fsa_fname);
	$fsa = keepACGTonly($fsa);
	my $numseqs = scalar(@$fsa);

	print STDERR "CMD: $0 " . join(' ', @ARGV) . "\n";
	print STDERR "\n";
	print STDERR "Number of sequences: $numseqs\n";
	print_marker_params(*STDERR, 'Left', $params->{leftmark});
	print_marker_params(*STDERR, 'Mid', $params->{midmark});
	print_marker_params(*STDERR, 'Right', $params->{rightmark});
	print STDERR "\n";
	printf STDERR ("Maximum number of mismatches: %d\n", $params->{max_mismatch});
	printf STDERR ("Offset for left primer: %d\n", $params->{offsetleft});
	printf STDERR ("Offset for right primer: %d\n", $params->{offsetright});
	print STDERR "\n";
	my $trimneg = $params->{trimneg};
	my $trimpos = $params->{trimpos};
	my $longlen_trimneg = $params->{longlen_trimneg};
	my $longlen_trimpos = $params->{longlen_trimpos};
	if(defined($trimneg) && defined($trimpos)) {
		print STDERR "trim: -$trimneg, +$trimpos\n";
	}
	if(defined($longlen_trimneg) && defined($longlen_trimpos)) {
		print STDERR "longlen-trim: -$longlen_trimneg, +$longlen_trimpos\n";
	}
	print STDERR "\n";

	if(length($params->{midmark}->{str}) == 0) {
		die "middle primer must be defined";
	}

	#sequence arr 
	$params->{leftmark}->{arr} = consensus2arr($params->{leftmark});
	$params->{rightmark}->{arr} = consensus2arr($params->{rightmark});
	$params->{midmark}->{arr} = consensus2arr($params->{midmark});

	my @new_fsa = ();
	my @longlen_filtered_fsa = ();
	foreach my $s (@$fsa) {
		#pivot is the position of mid-str
		my ($has_found, $pivot) = isolate_region_helper(
			*STDERR, 
			$s->{tag},
			$s->{seq}, 
			$params, 
		);
		my $seqlen = length($s->{seq});

		my $is_good = 1;
		my $shortseq = '';
		if(!$has_found) {
			$is_good = 0;
			printf STDERR ("Removed because of missing markers\n");
		}
		if($is_good) {
			if($pivot >= $trimneg && $pivot + $trimpos < $seqlen ) {
				$shortseq = substr($s->{seq}, $pivot-$trimneg, $trimneg + $trimpos);
			}
			else {
				$is_good = 0;
				printf STDERR ("Removed because of truncation. seqlen=$seqlen\n");
			}
		}

		#trim the original to a shorter length as well
		my $longseq = $s->{seq};
		if($is_good && defined($longlen_trimneg) && defined($longlen_trimpos)) {
			my $long_startpos = $pivot - $longlen_trimneg;
			my $long_len = $longlen_trimneg + $longlen_trimpos;
			if($long_startpos >= 0 && $long_startpos + $long_len < $seqlen) {
				$longseq = substr($s->{seq}, $long_startpos, $long_len);
			}
			else {
				$is_good = 0;
				printf STDERR ("Removed because longlen-trim. seqlen=" .length($s->{seq}). "\n");
			}
		}

		if($is_good) {
			my %isolated_hash = (
				tag => $s->{tag}, 
				seq => $shortseq,
			);
			my %longseq_hash = (
				tag => $s->{tag}, 
				seq => $longseq,
			);

			push(@new_fsa, \%isolated_hash);
			push(@longlen_filtered_fsa, \%longseq_hash); #keep the original length
		}
		else {
			printf STDERR ("Removed: %s\n", $s->{tag});
		}
		print STDERR "\n";
	}
	writeFa(*STDOUT, \@new_fsa);

	my $longlen_fname = $params->{longlen_fname};
	if(defined($longlen_fname)) {
		open(my $outfh, ">$longlen_fname") or die "Cannot open $longlen_fname for write";
		writeFa($outfh, \@longlen_filtered_fsa);
		close($outfh);
	}

	my @lenlst = ();
    foreach my $s (@$fsa) {
		push(@lenlst, length($s->{seq}));
	}
	@lenlst = sort {$a <=> $b} @lenlst;
	printf STDERR ("Min sequence length: %d\n", $lenlst[0]);
	printf STDERR ("5%%-quantile sequence length: %d\n", compute_empirical_quantile(\@lenlst, 0.05));
	printf STDERR ("Median sequence length: %d\n", compute_empirical_quantile(\@lenlst, 0.50));
	printf STDERR ("95%%-quantile sequence length: %d\n", compute_empirical_quantile(\@lenlst, 0.95));
	printf STDERR ("Max sequence length: %d\n", $lenlst[scalar(@lenlst)-1]);
}

sub print_marker_params
{
	my ($fptr, $tag, $marker) = @_;

	if($marker->{isValid}) {
		printf $fptr ("$tag k-mer/coord: %s (%d)\n", 
			$marker->{str}, $marker->{coord},
		);
	}
	else {
		printf $fptr ("$tag k-mer/coord: undefined\n");
	}
}

sub create_nullmarker
{
	#str, coord, arr, isValid
	my %hash = (
		str => undef, 
		coord => undef,
		arr => undef, 
		isValid => 0,
	);
	return \%hash;
}


sub isolate_region_helper
{
	my (
		$fptr,
		$seqname,
		$seq, 
		$params,
	) = @_;

	print STDERR "$seqname\n";

	$seq = uc($seq);
	$seq =~ s/[^ACGT]//g;
	#$seq =~ tr/ACGT/0123/; #not compatible with R and Y characters
	my @seqarr = split(//, $seq);

	my $leftlst = find_indices_of_matches($fptr, $params->{leftmark}{arr}, \@seqarr, $params->{max_mismatch});
	my $midlst = find_indices_of_matches($fptr, $params->{midmark}{arr}, \@seqarr, $params->{max_mismatch});
	my $rightlst = find_indices_of_matches($fptr, $params->{rightmark}{arr}, \@seqarr, $params->{max_mismatch});

	my @recordlst = ();
	foreach my $mid (@$midlst) {
		my $is_goodpivot = 1;
		my %record = (
			left => -1, 
			right => -1, 
			mid => $mid,
		);
		if($params->{leftmark}{isValid}) {
			my $found_left = 0;
			my $lit_diff = $params->{midmark}->{coord} - $params->{leftmark}->{coord};
			foreach my $left (@$leftlst) {
				if(
					$left < $mid 
					&& abs( ($mid-$left) - $lit_diff) <= $params->{offsetleft}
				) {
					$found_left = 1;
					$record{left} = $left;
					last;
				}
			}
			if(!$found_left) {
				$is_goodpivot = 0;
			}
		}
		if($is_goodpivot && $params->{rightmark}{isValid}) {
			my $found_right = 0;
			my $lit_diff = $params->{rightmark}->{coord} - $params->{midmark}->{coord};
			foreach my $right (@$rightlst) {
				if(
					$right > $mid 
					&& abs( ($right -$mid) - $lit_diff) <= $params->{offsetright}
				) {
					$found_right = 1;
					$record{right} = $right;
					last;
				}
			}
			if(!$found_right) {
				$is_goodpivot = 0;
			}
		}
		if($is_goodpivot) {
			push(@recordlst, \%record);
		}
	}

	#display everything on the recordlst
	foreach my $record (@recordlst) {
		printf STDERR ("left=%d  mid=%d  right=%d\n", 
			$record->{left}, $record->{mid}, $record->{right},
		);
	}
	print STDERR "Number of matches satisfying markers: " . scalar(@recordlst) . "\n";

	if(scalar(@recordlst) >= 1) {
		return (1, $recordlst[0]->{mid});
	}
	else {
		return (0, undef);
	}
}

sub compute_log_pwm
{
	my $pwm = shift;

	foreach my $v (@$pwm) {
		for(my $a = 0; $a < scalar(@$v); $a++) {
			$v->[$a] = log($v->[$a]);
		}
	}
}


sub add_pseudofreq
{
	my ($pwm, $ntfreq, $pseudoweight) = @_;
	foreach my $v (@$pwm) {
		for(my $a = 0; $a < scalar(@$v); $a++) {
			$v->[$a] = $v->[$a] * (1.0 - $pseudoweight) + $ntfreq->[$a] * $pseudoweight;
		}
	}
}

sub compute_ntfreq
{
	my ($seqarr) = @_;

	my @count = (0, 0, 0, 0);

	foreach my $a (@$seqarr) {
		$count[$a]++;
	}
	my $sum = sum(@count);

	for(my $a = 0; $a < scalar(@count); $a++) {
		$count[$a] /= $sum;
	}

	return \@count;
}



