#!/usr/bin/perl -w
use strict;

package PwmDistance;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(matrix_distance print_pwm calc_motif_distance1 calc_motif_distance2 normalize_pwm);

#Type I - minimum overlap width is 7. If the smaller motif is shorter than 7 bases, 
#         the minimum is 2 bases fewer than the smaller motif length.
#Type II - minimum overlap width is 6. If the smaller motif is shorter than 6 bases, 
#          all bases of the smaller motif are used.
#          Additional constraint that the average entropy of the learned motif must 
#          be at least 1 in the overlapping region

####################################################################################
# Type I 
####################################################################################

sub helper_calc_motif_distance1
{
	my ($pwm1, $pwm2, $normalization_mode) = @_;	


	my $dist1 = calc_motif_distance_forward_strand1_v2($pwm1, $pwm2, $normalization_mode);
	my $dist2 = calc_motif_distance_forward_strand1_v2($pwm1, pwm_revcompl($pwm2), $normalization_mode);

	#take min()
	if($dist1 < $dist2) {
		return $dist1;
	}
	else {
		return $dist2;
	}
}
sub calc_motif_distance1
{
	#normalization_mode: rms is the correct rms, buggy_mse is the version from the original Harbison suppl method section
	my ($pwm1, $pwm2, $normalization_mode);	
	if(scalar(@_) == 3) {
		($pwm1, $pwm2, $normalization_mode) = @_;	
	}
	else {
		($pwm1, $pwm2) = @_;	
		$normalization_mode = 'rms';
	}

	return helper_calc_motif_distance1($pwm1, $pwm2, $normalization_mode);
}

sub calc_motif_distance_forward_strand1_v2
{
	my ($pwm1, $pwm2, $normalization_mode) = @_;	

	#is_prob_matrix($pwm1);
	#is_prob_matrix($pwm2);

	my $smaller_len = (scalar(@$pwm1) < scalar(@$pwm2) ? scalar(@$pwm1) : scalar(@$pwm2));
	my $min_overlap_len = $smaller_len - 2;
	if($min_overlap_len > 7) {
		$min_overlap_len = 7;
	}

	my $min_dist = 999999.99;

	for(my $i = 0; $i < scalar(@$pwm1); $i++) {
		for(my $j = 0; $j < scalar(@$pwm2); $j++) {
			if($i == 0 || $j == 0) {
				my $align_len;
				if(scalar(@$pwm1) - $i < scalar(@$pwm2) - $j) {
					$align_len = scalar(@$pwm1) - $i;
				}
				else {
					$align_len = scalar(@$pwm2) - $j;
				}
				if($align_len >= $min_overlap_len) {
					my $dist = matrix_distance_with_offset($pwm1, $pwm2, $i, $j, $align_len, $normalization_mode);
					$min_dist = ($min_dist < $dist ? $min_dist : $dist); # min()
				}
			}
		}
	}
	return $min_dist;
}

sub calc_motif_distance_forward_strand1
{
	my ($pwm1, $pwm2, $normalization_mode) = @_;	

	is_prob_matrix($pwm1);
	is_prob_matrix($pwm2);

	#print_pwm($pwm1);
	#print_pwm($pwm2);
	
	my ($small, $large);
	if(scalar(@$pwm1) < scalar(@$pwm2)) {
		$small = $pwm1;
		$large = $pwm2;
	}
	else {
		$small = $pwm2;
		$large = $pwm1;
	}

	my $min_dist = 999999.99;

	my $min_overlap_len = scalar(@$small) - 2;
	if($min_overlap_len > 7) {
		$min_overlap_len = 7;
	}

	# Case 1: SMALL has variable starting position and LARGE has 0
	# XXXX
	#   YYYYYYYYYYYYY 
	# Here $i is 2 because the starting position of SMALL is 2. 
	# Starting position is always 0 for LARGE
	for(my $i = 0; $i < scalar(@$small); $i++) {
		my $remain_small_len = scalar(@$small) - $i;

		#Following line is not needed because it $align_len should always be set to $remain_small_len
		#Only have it for consistency with case 2
		my $align_len = ($remain_small_len < scalar(@$large) ? $remain_small_len : scalar(@$large)); #min()

		if($align_len >= $min_overlap_len) {
			my @align1 = ();
			my @align2 = ();
			for(my $k = 0; $k < $align_len; $k++) {
				push(@align1, $small->[$i + $k]);
				push(@align2, $large->[$k]);
			}
			my $dist = matrix_distance(\@align1, \@align2, $normalization_mode);

			#print "case1, pos $i, $dist\n";
			$min_dist = ($min_dist < $dist ? $min_dist : $dist); # min()
		}
	}

	# Case 2: LARGE has variable starting position and SMALL has 0
	for(my $i = 0; $i < scalar(@$large); $i++) {
		my $remain_large_len = scalar(@$large) - $i;
		my $align_len = ($remain_large_len < scalar(@$small) ? $remain_large_len : scalar(@$small)); #min()

		if($align_len >= $min_overlap_len) {
			my @align1 = ();
			my @align2 = ();
			for(my $k = 0; $k < $align_len; $k++) {
				push(@align1, $small->[$k]);
				push(@align2, $large->[$i + $k]);
			}
			my $dist = matrix_distance(\@align1, \@align2, $normalization_mode);
			#print "case2, pos $i, $dist\n";
			$min_dist = ($min_dist < $dist ? $min_dist : $dist); # min()
		}
	}
	return $min_dist;

}
####################################################################################
# Type II
####################################################################################
sub calc_motif_distance2
{
	my ($predict_pwm, $target_pwm, $normalization_mode);

	$normalization_mode = 'rms';
	my $is_modify = 1; #use our modified version by ignoring N chars
	($predict_pwm, $target_pwm, $is_modify, $normalization_mode) = @_;	

	return helper_calc_motif_distance2($predict_pwm, $target_pwm, $is_modify, $normalization_mode);
}



sub helper_calc_motif_distance2
{
	my ($predict_pwm, $target_pwm, $is_modify, $normalization_mode) = @_;

	my $dist1 = calc_motif_distance_forward_strand2($predict_pwm, $target_pwm, $is_modify, $normalization_mode);
	my $dist2 = calc_motif_distance_forward_strand2($predict_pwm, pwm_revcompl($target_pwm), $is_modify, $normalization_mode);

	#take min()
	if($dist1 < $dist2) {
		return $dist1;
	}
	else {
		return $dist2;
	}
}

sub calc_motif_distance_forward_strand2
{
	my ($predict, $target, $is_modify, $normalization_mode) = @_;	

	is_prob_matrix($predict);
	is_prob_matrix($target);

	my $min_dist = 999999.99;

	my $min_overlap_len = 6;
	if(scalar(@$predict) < 6 || scalar(@$target) < 6) {
		#take the min of the two motif length
		if(scalar(@$predict) < scalar(@$target)) {
			$min_overlap_len = scalar(@$predict);
		}
		else {
			$min_overlap_len = scalar(@$target);
		}
	}

	# Case 1: PREDICT has variable starting position and TARGET has 0
	# XXXX
	#   YYYYYYYYYYYYY 
	# Here $i is 2 because the starting position of PREDICT is 2. 
	# Starting position is always 0 for TARGET
	for(my $i = 0; $i < scalar(@$predict); $i++) {
		my $remain_len = scalar(@$predict) - $i;

		my $align_len = ($remain_len < scalar(@$target) ? $remain_len : scalar(@$target)); #min()

		if($align_len >= $min_overlap_len) {
			my @align_p = ();
			my @align_t = ();
			for(my $k = 0; $k < $align_len; $k++) {
				push(@align_p, $predict->[$i + $k]);
				push(@align_t, $target->[$k]);
			}

			#whether statisfy avg_entropy condition
			my $is_satisfy_avg_entropy_cond = 0;
			if($is_modify) {
				if(avg_entropy(remove_n_columns(\@align_p, \@align_t)) >= 1.0) {
					$is_satisfy_avg_entropy_cond = 1;
				}
			}
			else {
				if(avg_entropy((\@align_p)) >= 1.0) {
					$is_satisfy_avg_entropy_cond = 1;
				}
			}
			
			if($is_satisfy_avg_entropy_cond) {
				my $dist = matrix_distance(\@align_p, \@align_t, $normalization_mode);
				#print "case1, pos $i, $dist\n";
				$min_dist = ($min_dist < $dist ? $min_dist : $dist); # min()
			}
		}
	}

	# Case 2: TARGET has variable starting position and PREDICT has 0
	for(my $i = 0; $i < scalar(@$target); $i++) {
		my $remain_len = scalar(@$target) - $i;
		my $align_len = ($remain_len < scalar(@$predict) ? $remain_len : scalar(@$predict)); #min()

		if($align_len >= $min_overlap_len) {
			my @align_p = ();
			my @align_t = ();
			for(my $k = 0; $k < $align_len; $k++) {
				push(@align_p, $predict->[$k]);
				push(@align_t, $target->[$i + $k]);
			}
			#whether statisfy avg_entropy condition
			my $is_satisfy_avg_entropy_cond = 0;
			if($is_modify) {
				if(avg_entropy(remove_n_columns(\@align_p, \@align_t)) >= 1.0) {
					$is_satisfy_avg_entropy_cond = 1;
				}
			}
			else {
				if(avg_entropy((\@align_p)) >= 1.0) {
					$is_satisfy_avg_entropy_cond = 1;
				}
			}

			if($is_satisfy_avg_entropy_cond) {
				my $dist = matrix_distance(\@align_p, \@align_t, $normalization_mode);
				#print "case2, pos $i, $dist\n";
				$min_dist = ($min_dist < $dist ? $min_dist : $dist); # min()
			}
		}
	}
	return $min_dist;

}

##########################################################################################
# Helper functions
##########################################################################################

sub remove_n_columns {
	my ($align_predict, $align_target) = @_;

	if(scalar(@$align_predict) != scalar(@$align_target)) {
		die "Matrices are not aligned";
	}

	my @new_align_predict = ();
	for(my $i = 0; $i < scalar(@$align_predict); $i++) {
		my $is_N = 0; # all positions of target are 0.25
		foreach my $f (@{$align_target->[$i]}) {
			if( ($f - 0.25) < 0.0001) {
				$is_N = 1;
			}
		}
		if(!$is_N) {
			push(@new_align_predict, $align_predict->[$i]);
		}
	}
	return \@new_align_predict;

}

sub avg_entropy 
{
	my $pwm = shift;

	if(scalar(@$pwm) == 0) {
		return 2.0;
	}

	my $sum_entropy = 0.0;
	foreach my $v (@$pwm) {
		my $entropy = 2.0;
		foreach my $p (@$v) {
			if($p > 1e-30 ) {
				$entropy += $p * log($p) / log(2.0);
			}
		}
		$sum_entropy += $entropy;
	}
	my $avg = ($sum_entropy / scalar(@$pwm));
	return $avg;
}

sub is_prob_matrix
{
	my $pwm = shift;

	foreach my $vect (@$pwm) {
		my $sum = 0.0;
		if(scalar(@$vect) != 4) {
			die "Error: invalid prob matrix dimension";
		}
		foreach my $i (@$vect) {
			$sum += $i;
		}
		if(abs($sum - 1.0) > 0.000001) {
			die "Error: non-probability matrix found";
		}
	}
}

sub matrix_distance_with_offset
{
	my ($pwm1, $pwm2, $pos1, $pos2, $len, $normalization_mode) = @_;	

	my $result = 0.0;
	
	if($normalization_mode eq 'buggy_mse') {
		for(my $i = 0; $i < $len; $i++) {
			my $per_column_sum = 0.0;
			for(my $j = 0; $j < 4; $j++) {
				my $d = ($pwm1->[$pos1 + $i][$j] - $pwm2->[$pos2 + $i][$j]);
				$per_column_sum += $d * $d; #takes a square
			}
			$result += $per_column_sum;
		}
	}
	elsif($normalization_mode eq 'rms') {
		for(my $i = 0; $i < $len; $i++) {
			my $per_column_sum = 0.0;
			for(my $j = 0; $j < 4; $j++) {
				my $d = ($pwm1->[$pos1 + $i][$j] - $pwm2->[$pos2 + $i][$j]);
				$per_column_sum += $d * $d; #takes a square
			}
			$result += sqrt($per_column_sum);
		}
	}
	else {
		die "Invalid normalization_mode";
	}

	$result /= ($len * sqrt(2.0));
	return $result;
}

sub matrix_distance
{
	my ($pwm1, $pwm2, $normalization_mode) = @_;	

	#print_pwm($pwm1);
	#print_pwm($pwm2);

	if(scalar(@$pwm1) != scalar(@$pwm2)) {
		die "Error: different dimension in matrix_distance()";
	}
	my $span = scalar(@$pwm1);
	
	my $result = 0.0;
	for(my $i = 0; $i < $span; $i++) {
		my $per_column_sum = 0.0;
		for(my $j = 0; $j < 4; $j++) {
			$per_column_sum += ($pwm1->[$i][$j] - $pwm2->[$i][$j]) ** 2; #takes a square
		}
		if($normalization_mode eq 'buggy_mse') {
			$result += $per_column_sum;
		}
		elsif($normalization_mode eq 'rms') {
			$result += sqrt($per_column_sum);
		}
		else {
			die "Invalid normalization_mode";
		}
	}

	$result /= ($span * sqrt(2.0));
	return $result;
}

sub normalize_pwm
{
	my $pwm = shift;

	foreach my $vec (@$pwm) {
		my $sum = 0.0;
		for(my $i = 0; $i < scalar(@$vec); $i++) {
			$sum += $vec->[$i];
		}
		for(my $i = 0; $i < scalar(@$vec); $i++) {
			$vec->[$i] /= $sum;
		}
	}
}

sub pwm_revcompl 
{
	my $pwm = shift;

	#WARNING: need to make a copy first or the complement-step will screw up the original matrix
	my @new_pwm = ();
	for(my $i = 0; $i < scalar(@$pwm); $i++) {
		for(my $j = 0; $j < 4; $j++) {
			$new_pwm[$i][$j] = $pwm->[scalar(@$pwm)-$i-1][$j];
		}
	}

	#Complement-step
	#Swap cell 0 and 3 (A and T)
	#Swap cell 1 and 2 (C and G)
	my $tmp;
	for(my $i = 0; $i < scalar(@new_pwm); $i++) {
		#Swap cell 0 and 3 (A and T)
		$tmp = $new_pwm[$i][0];
		$new_pwm[$i][0] = $new_pwm[$i][3];
		$new_pwm[$i][3] = $tmp;

		#Swap cell 1 and 2 (C and G)
		$tmp = $new_pwm[$i][1];
		$new_pwm[$i][1] = $new_pwm[$i][2];
		$new_pwm[$i][2] = $tmp;
	}
	return \@new_pwm;
}

sub print_pwm
{
    my ($fh, $mat) = @_;

    for(my $i = 0; $i < scalar(@$mat); $i++) {
        printf $fh ("[%3d]", $i);
        for(my $j = 0; $j < scalar(@{$mat->[$i]}); $j++) {
            printf $fh (" %8.3f", $mat->[$i][$j]);
        }
        print $fh "\n";
    }
}


1;
__END__
