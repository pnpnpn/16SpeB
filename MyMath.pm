#!/usr/bin/perl -w
use strict;


package MyMath;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
	nchoosek
	round
	compute_empirical_quantile
	compute_avg_stdev
	compute_rank
	change_base
	sign
	log10
	log2
	mat_transpose
);

use POSIX qw(ceil floor);
use List::Util qw(sum min max);

sub mat_transpose
{
	my $mat = shift;

	my $dim2 = scalar(@{$mat->[0]});
	for(my $i = 0; $i < scalar(@$mat); $i++) {
		if(scalar(@{$mat->[$i]}) != $dim2) {
			die "inconsistent matrix dimension found in mat_transpose()";
		}
	}
	my @new_mat = ();

	for(my $i = 0; $i < scalar(@$mat); $i++) {
		for(my $j = 0; $j < scalar(@{$mat->[$i]}); $j++) {
			$new_mat[$j][$i] = $mat->[$i][$j];
		}
	}

	return \@new_mat;
}

sub sign {
	my ($n) = @_;
	if($n >= 0) {
		return +1;
	}
	else {
		return -1;
	}
}


sub compute_rank {
	my ($v, $lst) = @_;

	my @sorted = sort {$a <=> $b} @$lst;
	for(my $i = 0; $i < scalar(@sorted); $i++) {
		if($v < $sorted[$i]) {
			return $i;
		}
	}
	return scalar(@sorted);
}

sub compute_avg_stdev {
	my $lst = shift;

	my $avg= sum(@$lst) / scalar(@$lst);

	my $std= 0;
	foreach my $l (@$lst) {
		my $d = $l - $avg;
		$std+= $d * $d;
	}
	$std = sqrt($std / (scalar(@$lst) - 1));
	return ($avg, $std);
}


sub round
{
	my $number = shift;
	return floor($number + .5);
}

sub log10 {
	my $x = shift;
	return (log($x) / log(10.0));
}
sub log2 {
	my $x = shift;
	return (log($x) / log(2.0));
}

sub change_base
{
    my ($decimal, $base, $len) = @_;

	#1 is converted to 0..001
	#2 is converted to 0..010

    if($decimal >= $base ** $len) {
        die "Invalid usage of change_base: len=$len, decimal=$decimal, base=$base";
    }

    my @arr = ();
    for(my $i = 0; $i < $len; $i++) {
        $arr[$i] = 0; #initialize
    }

    my $i = 0;
    while($decimal != 0) {
        $arr[($len-1) -$i] = $decimal % $base;
        $decimal = floor($decimal / $base);
        $i++;
    }
    return \@arr;
}

sub compute_empirical_quantile
{
    my ($sorted, $q) = @_;

    my $numsamples = scalar(@$sorted);
    if ($numsamples == 0) {
        die "Invalid numsamples=0\n";
    }

    my $integral = floor(($numsamples-1) * $q);
    my $fraction = ($numsamples-1) * $q - $integral;

    if($integral < 0 || $integral >= $numsamples) {
        die "Invalid integral $integral\n";
    }
    if($fraction < 0.0 || $fraction >= 1.0 ) {
        die "Invalid fractional $fraction\n";
    }

    if($integral + 1 >= $numsamples) {
        return $sorted->[$integral];
    }
    else {
        #interpolate
        return $sorted->[$integral] + $fraction * ($sorted->[$integral+1] - $sorted->[$integral]);
    }
}



sub nchoosek
{
	my ($n, $k) = @_;

	my $sum = 0.0;

	#compute n! / (n-k)!
	for(my $i = $n-$k+1; $i <= $n; $i++) {
		$sum += log($i);
	}

	#compute k!
	for(my $i = 1; $i <= $k; $i++) {
		$sum -= log($i);
	}

	return round(exp($sum));
}


1;

__END__
