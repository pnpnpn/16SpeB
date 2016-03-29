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
use ConsensusSeq;
use MyRandom;
use MyMath;
use MatrixPermutation;
use PwmDistance;
use POSIX qw(ceil floor);
use List::Util qw(sum min max);

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
usage: <program> --fsa=<FSA> --left-str=<STRING> --right=str=<STRING> --left-range=<INT,INT> --right-range=<INT,INT>

Isolates a region where endpoints are best matches of 2 user-specified k-mers.

OPTIONS:
left-str                 left k-mer
right-str                right k-mer
pw=<FLT>                 pseudoweight (default: 0.1)
left-llr=<FLT>           LLR cutoff for dropping a sequence
right-llr=<FLT>          LLR cutoff for dropping a sequence
len-left=<INT>           count len from left and then truncate
len-right=<INT>          count len from right and then truncate

longlen-fsa=<FILE>        writes to file the filtered seqs with of longer length
longlen-rtrim-pos=<INT>   trim in the positive direction based on the right-str
longlen-rtrim-neg=<INT>   trim in the negative direction based on the right-str

example: r16s_scripts/isolate_region.pl --fsa=Lactobacillus_seqpair.fsa --left-str=AGAGTTTGATCCTGGCTCAG --right-str=ACTCCTACGGGAGGCAGCAY --left-range=0,0 --right-range=250,600 --len-right=270 --longlen-rtrim-pos=1000 --longlen-rtrim-neg=270 --longlen-fsa=filter.fsa

Created on 6/03/2010

USAGE
		exit(1);
	}

	my $fsa_fname = undef;
	my $pseudoweight = 0.1;
	my $left_mer = undef;
	my $right_mer = undef;
	my $left1 = undef;
	my $left2= undef;
	my $right1= undef;
	my $right2= undef;
	my $left_llr_cutoff= undef;
	my $right_llr_cutoff= undef;
	my $truncate_len = undef;
	my $truncate_ref = undef;
	my $longlen_fname= undef;
	my $longlen_rtrim_pos= undef;
	my $longlen_rtrim_neg= undef;
	foreach my $a (@ARGV) {
		if( $a =~ /^--left-str=(\S+)/) {
			$left_mer = $1;
		}
		elsif( $a =~ /^--right-str=(\S+)/) {
			$right_mer = $1;
		}
		elsif( $a =~ /^--left-range=(\d+),(\d+)/) {
			$left1= $1;
			$left2= $2;
		}
		elsif( $a =~ /^--right-range=(\d+),(\d+)/) {
			$right1= $1;
			$right2= $2;
		}
		elsif( $a =~ /^--fsa=(\S+)/) {
			$fsa_fname=$1;
		}
		elsif( $a =~ /^--pw=([\d\.]+)/) {
			$pseudoweight = $1;
		}
		elsif( $a =~ /^--left-llr=([\-\d\.]+)/) {
			$left_llr_cutoff= $1;
		}
		elsif( $a =~ /^--right-llr=([\-\d\.]+)/) {
			$right_llr_cutoff= $1;
		}
		elsif( $a =~ /^--len-left=(\d+)/) {
			$truncate_len = $1;
			$truncate_ref = 'left';
		}
		elsif( $a =~ /^--len-right=(\d+)/) {
			$truncate_len = $1;
			$truncate_ref = 'right';
		}
		elsif( $a =~ /^--longlen-fsa=(\S+)/) {
			$longlen_fname = $1;
		}
		elsif( $a =~ /^--longlen-rtrim-neg=(\d+)/) {
			$longlen_rtrim_neg= $1;
		}
		elsif( $a =~ /^--longlen-rtrim-pos=(\d+)/) {
			$longlen_rtrim_pos= $1;
		}
		else {
			die "Unrecognized parameter: $a\n";
		}
	}
	if(defined($longlen_rtrim_neg) xor defined($longlen_rtrim_pos)) {
		die "Both longlen-trim pos/neg must be defined";
	}
	if($longlen_rtrim_neg < 0 || $longlen_rtrim_pos < 0) {
		die "trim pos/neg must be positive";
	}
	
	#if a new random seed is read
	#if(defined($user_randseed)) {
	#	$RAND_SEED = $user_randseed;
	#}
	$RAND_MT_GEN->srand($RAND_SEED);
	#print STDERR "Random seed: " . $RAND_SEED . "\n";

	my $fsa = parseFa($fsa_fname);
	my $numseqs = scalar(@$fsa);

	print STDERR "CMD: $0 " . join(' ', @ARGV) . "\n";
	print STDERR "\n";
	print STDERR "Number of sequences: $numseqs\n";
	print STDERR "Left k-mer: $left_mer\n";
	print STDERR "Right k-mer: $right_mer\n";
	print STDERR "Left range: $left1,$left2\n";
	print STDERR "Right range: $right1,$right2\n";
	print STDERR "Pseudo-weight: $pseudoweight\n";
	if(defined($left_llr_cutoff)) {
		print STDERR "Left LLR-threshold: $left_llr_cutoff\n";
	}
	if(defined($right_llr_cutoff)) {
		print STDERR "Right LLR-threshold: $right_llr_cutoff\n";
	}
	print STDERR "\n";
	if(defined($truncate_ref)){
		print STDERR "Trucate point of reference: $truncate_ref\n";
		print STDERR "Isolated trucate length: $truncate_len\n";
	}
	if(defined($longlen_rtrim_neg) && defined($longlen_rtrim_pos)) {
		print STDERR "longlen-rtrim: -$longlen_rtrim_neg, +$longlen_rtrim_pos\n";
	}
	print STDERR "\n";


	#construct matrix from k-mer
	my $left_pwm = consensus2pwm($left_mer, 0.0);
	my $right_pwm = consensus2pwm($right_mer, 0.0);

	print STDERR "Left PWM\n";
	print_pwm(*STDERR, $left_pwm);
	print STDERR "\n";
	print STDERR "Right PWM\n";
	print_pwm(*STDERR, $right_pwm);
	print STDERR "\n";

	my @new_fsa = ();
	my @longlen_filtered_fsa = ();
	foreach my $s (@$fsa) {
		my ($newseq, $lllr, $rllr, $lpivot, $rpivot) = isolate_region_helper(
			*STDERR, 
			$s->{tag},
			$s->{seq}, 
			$left_pwm, 
			$right_pwm, 
			$left1,
			$left2,
			$right1,
			$right2,
			$pseudoweight,
		);

		my $is_good = 1;
		if(defined($left_llr_cutoff) && $lllr < $left_llr_cutoff) {
			$is_good = 0;
			printf STDERR ("Removed because of primer not found\n");
		}
		if(defined($right_llr_cutoff) && $rllr < $right_llr_cutoff) {
			$is_good = 0;
			printf STDERR ("Removed because of primer not found\n");
		}
		if($is_good && defined($truncate_ref)) {
			if(length($newseq) < $truncate_len) {
				$is_good = 0;
				printf STDERR ("Removed because of truncation\n");
			}
			else {
				my $curlen = length($newseq);
				$newseq = substr($newseq, $curlen - $truncate_len, $truncate_len);
			}
		}

		#trim the original to a shorter length as well
		my $long_seq = $s->{seq};
		if($is_good && defined($longlen_rtrim_neg) && defined($longlen_rtrim_pos)) {
			my $long_startpos = $rpivot - $longlen_rtrim_neg;
			my $long_len = $longlen_rtrim_neg + $longlen_rtrim_pos;
			if($long_startpos >= 0 && $long_startpos + $long_len < length($s->{seq})) {
				$long_seq = substr($s->{seq}, $long_startpos, $long_len);
			}
			else {
				$is_good = 0;
				printf STDERR ("Removed because longlen-trim. seqlen=" .length($s->{seq}). "\n");
			}
		}

		if($is_good) {
			my %isolated_hash = (
				tag => $s->{tag}, 
				seq => $newseq,
			);
			my %longseq_hash = (
				tag => $s->{tag}, 
				seq => $long_seq,
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


sub isolate_region_helper
{
	my (
		$fptr,
		$seqname,
		$seq, 
		$orig_left_pwm, 
		$orig_right_pwm, 
		$left1,
		$left2,
		$right1,
		$right2,
		$pseudoweight,
	) = @_;

	print STDERR "$seqname\n";

	my $leftpwm = copy_matrix($orig_left_pwm);
	my $rightpwm = copy_matrix($orig_right_pwm);

	$seq = uc($seq);
	$seq =~ s/[^ACGT]//g;
	$seq =~ tr/ACGT/0123/;
	my @seqarr = split(//, $seq);
	my $ntfreq = compute_ntfreq(\@seqarr);

	add_pseudofreq($leftpwm, $ntfreq, $pseudoweight);
	add_pseudofreq($rightpwm, $ntfreq, $pseudoweight);

	#display
	print $fptr "Background freq: ";
	foreach my $freq (@$ntfreq) {
		printf $fptr ("%.2lf ", $freq);
	}
	print $fptr "\n";
	#print_pwm(*STDERR, $leftpwm);
	#print_pwm(*STDERR, $rightpwm);

	compute_log_pwm($leftpwm);
	compute_log_pwm($rightpwm);

	my ($leftind, $left_llr) = find_ind_of_best_match($fptr, \@seqarr, $leftpwm, $left1, $left2);
	my ($rightind, $right_llr)  = find_ind_of_best_match($fptr, \@seqarr, $rightpwm, $right1, $right2);

	print $fptr "Best positions left=$leftind right=$rightind\n";

	if($leftind >= $rightind) {
		die "inconsistent left/right indices: $leftind and $rightind";
	}

	my $newseq = join('', @seqarr[$leftind .. $rightind -1]);
	$newseq =~ tr/0123/ACGT/;

	return ($newseq, $left_llr, $right_llr, $leftind, $rightind);
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

sub find_ind_of_best_match
{
	my ($fptr, $seqarr, $log_pwm, $range1, $range2) = @_;

	my $seqlen = scalar(@$seqarr);
	my $width = scalar(@$log_pwm);

	my $leftp = 0;
	my $rightp = $seqlen - $width + 1; #exclusive

	$leftp = max($leftp, $range1);
	$rightp= min($rightp, $range2);

	if($range1 >= $seqlen - $width) {
		die "range1 is too small: $range1, $seqlen, $width";
	}
	if($leftp > $rightp) {
		die "inconsistent left/right ranges: $leftp and $rightp";
	}

	my @scorearr = ();
	for(my $pos = $leftp; $pos < $rightp; $pos++) {
		$scorearr[$pos] = 0.0;
		for(my $i = 0; $i < $width; $i++) {
			$scorearr[$pos] += $log_pwm->[$i][ $seqarr->[$i+$pos] ];
		}
	}

	my $bestpos = $leftp;
	my $bestscore = $scorearr[$leftp];
	for(my $pos = $leftp; $pos < $rightp; $pos++) {
		if($bestscore < $scorearr[$pos]) {
			$bestscore = $scorearr[$pos];
			$bestpos = $pos;
		}
	}

	#print the best match k-mer
	my $str = join('', @$seqarr[$bestpos .. $bestpos+$width-1]);
	$str =~ tr/0123/ACGT/;
	if(defined($bestscore)) {
		print $fptr "Best match: $str\n";
		printf $fptr ("Best score: %.5lf\n", $bestscore);
	}

	return ($bestpos, $bestscore);
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



