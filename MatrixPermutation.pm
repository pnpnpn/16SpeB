#!/usr/bin/perl -w
use strict;
use File::Copy;
#use Math::Random;
use File::Basename;

package MatrixPermutation;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(copy_matrix);

sub permutate_matrix
{
	my $mat = shift;
	my $copy = copy_matrix($mat);
	return random_permutate_array($copy);
}

sub copy_matrix
{
	my $src = shift;

	my @dest = ();

	for(my $i = 0; $i < scalar(@$src); $i++) {
		for(my $j = 0; $j < scalar(@{$src->[$i]}); $j++) {
			$dest[$i][$j] = $src->[$i][$j];
		}
	}
	return \@dest;
}

#this is done by swapping
sub random_permutate_array
{
	my $arr = shift;
	
	my $len = scalar(@$arr);
	for(my $i = 0; $i < $len; $i++) {
		#my $rand = int(Math::Random::random_uniform() * ($len - $i));
		my $rand = int(rand() * ($len - $i));
		if($rand >= $len - $i) {
			die "Random number is too large: $rand at $i\n";
		}
		my $j = $rand + $i; #new index
		my $temp = $arr->[$i];
		$arr->[$i] = $arr->[$j];
		$arr->[$j] = $temp;
	}
}

1;
__END__

