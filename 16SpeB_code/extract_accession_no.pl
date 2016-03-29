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
usage: <program> <FSA-LIST>

Print out accession numbers for each file

Created on 4/4/2011

USAGE
		exit(1);
	}

	my $silva_count = 0;
	foreach my $fname (sort @ARGV) {
		if(! -f $fname) {
			die "File $fname does not exist";
		}
		my $fsa = parseFa($fname);
		my $numseqs = scalar(@$fsa);

		my @accno_lst = ();

		foreach my $s (@$fsa) {
			my $tag = $s->{tag};
			if($tag =~ /^([A-Za-z]+[\w\.]+)/) {
				push(@accno_lst, $1);
				$silva_count++;
			}
			elsif($tag =~ /^\d+ ([A-Za-z]+[\w\.]+)/) {
				push(@accno_lst, $1);
			}
			else {
				die "cannot parse tag: $tag";
			}
		}
		if(scalar(@accno_lst) != $numseqs) {
			die "inconsistent number of accession numbers";
		}

		#duplicate detection
		my %hash = ();
		foreach my $id (@accno_lst) {
			if(exists($hash{$id})) {
				print STDERR "duplicate id found: $id\n";
			}
			else {
				$hash{$id} = 1;
			}
		}

		printf ("Filename: $fname (%d sequences)\n", $numseqs);
		print join("\n", sort @accno_lst);
		print "\n================================================================\n\n";
	}
	print STDERR "Number of SILVA sequences: $silva_count\n";
}
