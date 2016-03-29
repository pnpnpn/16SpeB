#!/usr/bin/perl -w
use strict;
use diagnostics;

package R16sHelper;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
bin_strands_by_species
print_species_bins
find_indices_of_matches
consensus2arr
filter_duplicates_within_species_bins
);

use Benchmark qw(:all);
use File::Basename;
use File::Temp;
use FastaIO;
use FastaFilter;
use MyRandom;
use POSIX qw(ceil floor);

my $TMP_DIR = '/tmp/';

my $GAPCHAR = 4;
my $NUMALPHAS = 4;


sub consensus2arr
{
	my $marker= shift;
	my $str = $marker->{str};
	if(!$marker->{isValid}) {
		return undef;
	}
	if($str =~ /([^ACGTRYacgtry])/) {
		die "non-ACGTRY character found in $str";
	}
	$str = uc($str);
	#$str =~ tr/ACGT/0123/; #not compatible with R and Y characters
	my @arr = split(//, $str);
	return \@arr;
}

sub find_indices_of_matches
{
    my (
        $fptr,
        $kmer,
        $seqarr,
        $max_mismatch,
    ) = @_;

	if(!defined($kmer) || scalar(@$kmer) == 0) {
		return undef;
	}

    my @lst = ();
    for(my $i = 0; $i < scalar(@$seqarr) - scalar(@$kmer) + 1; $i++) {
        my $misscount = 0;
        for(my $j = 0; $j < scalar(@$kmer); $j++) {
            my $kmer_char = $kmer->[$j];
			my $seq_char = $seqarr->[$i + $j];
			my $ismatch = 0;
            if($kmer_char eq $seq_char) {
				$ismatch = 1;
            }
			elsif($kmer_char eq 'R') {
				if($seq_char =~ /[AG]/) { #purine is A or G
					$ismatch = 1;
				}
            }
			elsif($kmer_char eq 'Y') {
				if($seq_char =~ /[CT]/) { #pyrimidine is C or T
					$ismatch = 1;
				}
            }

			#bad match
			if(!$ismatch) {
				$misscount++;
			}
        }
        if($misscount <=  $max_mismatch) {
            push(@lst, $i);
        }
    }
    return \@lst;
}

sub bin_strands_by_species
{
	#fptr is only for outputing "stderr" output. 

	my ($fsa) = @_;
	my %species_bins = ();

	my $genus0 = undef;
	foreach my $s (@$fsa) {
		$s->{tag} =~ s/\"/ /g;
		my $tag = $s->{tag};
		$tag =~ s/ subsp\. / /g;
		$tag .= ' ';
		if(
			$tag =~ /^S\d+\s+(\w+) (\w+)[\s\;]/
			|| $tag =~ /^silva\|\S+\s+(\w+) (\w+)[\s\;]/
			|| $tag =~ /^\d+ \S+\s+(\w+) (\w+)[\s\;]/
		){
			my $genus = $1;
			my $species = $2;
			if(!exists($species_bins{$species})) {
				$species_bins{$species} = [];
			}
			push(@{$species_bins{$species}}, $s);

			#sanity check
			if(!defined($genus0)) {
				$genus0 = $genus;
			}
			if($genus0 ne $genus) {
				die "Unequal genus found $genus0 != $genus";
			}
		}
		else {
			print STDERR "WARNING: Cannot parse tag: " . $s->{tag} . "\n";
			die;
		}
	}
	return \%species_bins;
}

sub filter_duplicates_within_species_bins
{
	my ($species_bins) = @_;

	my %newbin = ();
	foreach my $species (keys %$species_bins) {
		$newbin{$species} = []; #must have at least one element
		my $lst = $species_bins->{$species};
		my %seen= ();
		foreach my $s (@$lst) {
			$s->{seq} = uc($s->{seq});
			$s->{seq} =~ s/[^ACGT]//g;
			if(!exists($seen{$s->{seq}})) {
				$seen{$s->{seq}} = 1;
				push(@{$newbin{$species}}, $s);

				#DEBUG
				if($species =~/oeni/) {
					print STDERR "darn! " . scalar(@$lst) . "\n";
					die;
				}

			}
			else {
				print STDERR "Strand removed because intra-species duplication " . $s->{tag} . "\n";
			}
		}
	}
	return \%newbin;
}

sub print_species_bins
{
	my ($fptr, $species_bins) = @_;
	my @species_names = sort keys %$species_bins;
	foreach my $name (@species_names) {
		my $num_strands = scalar(@{$species_bins->{$name}});
		my @taglst = ();
		foreach my $s (@{$species_bins->{$name}}) {
			my $tag = $s->{tag};
			$tag =~ s/\s/\_/g;
			push(@taglst, $tag);
		}
		printf $fptr ("%3d strands found for species $name: %s\n",
			$num_strands,
			join(" ", @taglst),
		);
	}
}

1;
__END__



