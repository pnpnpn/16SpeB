#!/usr/bin/perl -w
use strict;

package FastaFilter;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
	filter_by_genus
	filter_out_keywords
);


sub filter_by_genus
{
	my ($fptr, $fsa, $genus) = @_;

	my @new_fsa = ();
	foreach my $s (@$fsa) {
		if($s->{tag} =~ /$genus/i) {
			push(@new_fsa, $s);
		}
		else {
			print $fptr ("Removed (incorrect genus): ". $s->{tag} . "\n");
		}
	}
	print $fptr "\n";
	return \@new_fsa;
}

sub filter_out_keywords
{
	my ($fptr, $fsa, $keywords) = @_;
	my @new_fsa = ();
	foreach my $s (@$fsa) {
		my $found = 0;
		foreach my $word (@$keywords) {
			my $tag = $s->{tag} . " "; #to pick up words at the end of a line
			if($tag =~ /$word/i) {
				$found = 1;
			}
		}
		if(!$found) {
			push(@new_fsa, $s);
		}
		else {
			print $fptr ("Removed (bad keyword): ". $s->{tag} . "\n");
		}
	}
	print $fptr "\n";
	return \@new_fsa;
}



1;

__END__
