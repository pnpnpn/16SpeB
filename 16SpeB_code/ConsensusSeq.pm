#!/usr/bin/perl -w
use strict;

package ConsensusSeq;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(consensus2pwm pwm2consensus pwm2kmers_hack revcompl kmers2pwm nt2int int2nt);

use FastaIO;


#R	Purine (A or G)
#Y	Pyrimidine (C, T/U)
#M	C or A
#K	T/U or G
#W	T/U or A
#S	C or G

sub print_pwm
{
	my ($mat) = @_;

	for(my $i = 0; $i < scalar(@$mat); $i++) {
		printf ("%3d", $i+1);
		for(my $j = 0; $j < scalar(@{$mat->[$i]}); $j++) {
			printf (" %8.3f", $mat->[$i][$j]);
		}
		print "\n";
	}
}

sub consensus2pwm 
{
	my $consensus = undef;
	my $eps = undef;
	if(scalar(@_) == 1) {
		$consensus = shift;
		$eps = 0.012;
	}
	else {
		($consensus,$eps) = @_;
	}

	if(lc($consensus) =~ /[^acgtrymkwsn]/) {
		die "Invalid characters found at ($consensus)\n";
	}

	my $num_nts = 4;

	my @pwm = ();
	for(my $i = 0; $i < length($consensus); $i++) {
		for(my $j = 0; $j < $num_nts; $j++) {
			$pwm[$i][$j] = 0.0;
		}
	}

	my @arr = split(//, $consensus);
	for(my $i = 0; $i < scalar(@arr); $i++) {
		my $c = lc($arr[$i]);
		if($c =~ /[acgt]/i) {
			for(my $j = 0; $j < $num_nts; $j++) {
				if($j == nt2int($c)) {
					#$pwm[$i][$j] = 0.964;
					$pwm[$i][$j] = 1.0 - $eps;
				}
				else {
					$pwm[$i][$j] = $eps;
				}
			}
		}
		elsif($c eq 'n') {
			for(my $j = 0; $j < $num_nts; $j++) {
				$pwm[$i][$j] = 0.25;
			}
		}
		elsif($c =~ /[rymkws]/i) {
			for(my $j = 0; $j < $num_nts; $j++) {
				#$pwm[$i][$j] = 0.012;
				$pwm[$i][$j] = $eps;
			}
			if($c eq 'r') {
				$pwm[$i][nt2int('a')] = 0.5 - $eps;
				$pwm[$i][nt2int('g')] = 0.5 - $eps;
			}
			elsif($c eq 'y') {
				#$pwm[$i][nt2int('c')] = 0.488;
				#$pwm[$i][nt2int('t')] = 0.488;
				$pwm[$i][nt2int('c')] = 0.5 - $eps;
				$pwm[$i][nt2int('t')] = 0.5 - $eps;
			}
			elsif($c eq 'm') {
				$pwm[$i][nt2int('a')] = 0.5 - $eps;
				$pwm[$i][nt2int('c')] = 0.5 - $eps;
			}
			elsif($c eq 'k') {
				$pwm[$i][nt2int('g')] = 0.5 - $eps;
				$pwm[$i][nt2int('t')] = 0.5 - $eps;
			}
			elsif($c eq 'w') {
				$pwm[$i][nt2int('a')] = 0.5 - $eps;
				$pwm[$i][nt2int('t')] = 0.5 - $eps;
			}
			elsif($c eq 's') {
				$pwm[$i][nt2int('c')] = 0.5 - $eps;
				$pwm[$i][nt2int('g')] = 0.5 - $eps;
			}
		}
		else {
			die "Invalid characters found at ($consensus)\n";
		}
	}
	return \@pwm;
}

sub pwm2kmers_hack
{
	my ($pwm, $num_samples) = @_;
	#Generate kmers from a given PWM -- it is a hack because it simply samples from a PWM

	my %int2nt = (
		0 => 'A',
		1 => 'C',
		2 => 'G',
		3 => 'T',
	);

	my @kmers = ();
	for(my $i = 0; $i < $num_samples; $i++) {
		my @nts = ();
		for(my $m = 0; $m < scalar(@$pwm); $m++) {
			my $sum = 0.0; 
			for(my $a = 0; $a < scalar(@$pwm); $a++) {
				$sum += $pwm->[$m][$a];
			}
			my $rnd = rand() * $sum;

			my $cumsum = 0.0;
			my $nt = $int2nt{3}; #default
			for(my $a = 0; $a < scalar(@{$pwm->[$m]}); $a++) {
				$cumsum += $pwm->[$m][$a];
				if($rnd < $cumsum) {
					$nt = $int2nt{$a};
					last;
				}
			}
			$nts[$m] = $nt;
		}
		push( @kmers, join('', @nts) );
	}
	return \@kmers;
}

sub kmers2pwm
{
	my ($seqs) = @_;

	#ensure all lengths are the same
	my $len = length($seqs->[0]);
	foreach my $s (@$seqs) {
		if($len != length($s)) {
			die "Inconsistent length. All lengths must be equal: $s does not have length $len.";
		}
		if($s =~ /([^ACGTacgt]+)/) {
			die "Error: invalid character $1 in $s";
		}
	}

	#create matrix with dimension (motif-len, 4)
	my @mat = ();
	for(my $i = 0; $i < $len; $i++) {
		my $vec = [0, 0, 0, 0];
		push(@mat, $vec);
	}

	#tally
	foreach my $s (@$seqs) {
		$s = nt2int($s);
		my @sl = split(//, $s);
		for(my $i = 0; $i < $len; $i++) {
			$mat[$i][$sl[$i]]++;
		}
	}

	#normalize
	for(my $i = 0; $i < $len; $i++) {
		my $sum = 0; 
		for(my $j = 0; $j < scalar(@{$mat[$i]}); $j++) {
			$sum += $mat[$i][$j];
		}
		for(my $j = 0; $j < scalar(@{$mat[$i]}); $j++) {
			$mat[$i][$j] /= $sum;
		}
	}
	
	return \@mat;
}

sub pwm2consensus
{
	my ($freqmat) = @_;

	my $consensus = '';
	for(my $i = 0; $i < scalar(@$freqmat); $i++) {
		#"strong" is for the non-ambiguous nts (A,C,G or T)
		my @strong = ();
		#"weak" is for the ambiguous code (R,Y,M,K,W,S)
		my @weak = ();
		for(my $j = 0; $j < scalar(@{$freqmat->[$i]}); $j++) {
			if($freqmat->[$i][$j] > 0.5000001) {
				push(@strong, $j); #this should be unique
			}
			if($freqmat->[$i][$j] > 0.4) {
				push(@weak, $j); 
			}
		}


		#sanity check
		if(scalar(@weak) > 2) {
				die "Error: there is more than 2 elements in \@weak";
		}
		if(scalar(@strong) > 1) {
			die "Error: there is more than one element in \@strong. It should be unique.";
		}

		my $code;
		if(scalar(@weak) == 2) {
			$code = get_ambiguous_code(\@weak);
		}
		elsif(scalar(@strong) == 1) {
			$code = int2nt($strong[0]);
		}
		else {
			$code = 'n';
		}

		$consensus .= $code;
	}
	return $consensus;
}

sub revcompl_old
{
	my $nt_seq = shift;
	$nt_seq =~ tr/ACGTRYKMacgtrykm/TGCAYRMKtgcayrmk/;
	my $rev = reverse $nt_seq;
	return $rev;
}

sub get_ambiguous_code
{
	my ($weak) = @_;

	if(scalar(@$weak) != 2) {
		die "Error: nucleotide list must have exactly 2 elements.";
	}
	my ($nt1, $nt2);
	if($weak->[0] < $weak->[1]) {
		$nt1 = $weak->[0];
		$nt2 = $weak->[1];
	}
	else {
		$nt1 = $weak->[1];
		$nt2 = $weak->[0];
	}

	if($nt1 == nt2int('A') && $nt2 == nt2int('G')) {
		return 'r'; #purine A/G
	}
	elsif($nt1 == nt2int('C') && $nt2 == nt2int('T')) {
		return 'y'; #pyrimidine
	}
	elsif($nt1 == nt2int('G') && $nt2 == nt2int('T')) {
		return 'k'; #keto
	}
	elsif($nt1 == nt2int('A') && $nt2 == nt2int('C')) {
		return 'm'; #amino
	}
	elsif($nt1 == nt2int('C') && $nt2 == nt2int('G')) {
		return 's'; #strong
	}
	elsif($nt1 == nt2int('A') && $nt2 == nt2int('T')) {
		return 'w'; #weak
	}
	else {
		die "Invalid case in get_ambiguous_code -- nt1: $nt1, nt2: $nt2\n";
	}
}

sub nt2int
{
	my $str = shift;
	if($str =~ /[^ACGTacgt]/) {
		die "Invalid characters in $str at nt2int()";
	}
	$str =~ tr/ACGTacgt/01230123/;
	return $str;
}

sub int2nt
{
	my $str = shift;
	if($str =~ /[^0123]/) {
		die "Invalid characters in $str at int2nt()";
	}
	$str =~ tr/0123/ACGT/;
	return $str;
}

1;
__END__

