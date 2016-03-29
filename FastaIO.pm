#!/usr/bin/perl -w
use strict;
#use lib "$ENV{HOME}/perl_lib/usr/lib/perl5/site_perl/5.8.5/i386-linux-thread-multi/";

package FastaIO;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
	get_fasta_tmpdir
	revcompl 
	open_genome_file 
	validate_coord 
	open_genome_file_filter_coord 
	open_genome_file_keep_N
	parseFa 
	writeFa 
	isACGT
	keepACGTonly
);

use File::Copy;
use File::Temp qw(tempfile tempdir);
use File::Basename;


sub revcompl
{
	my $seq = shift;
	$seq = reverse($seq);
	#$seq =~ tr/ACGTacgt/TGCAtgca/;
    $seq =~ tr/ACGTRYKMacgtrykm/TGCAYRMKtgcayrmk/;
	return $seq;
}

sub get_fasta_tmpdir
{
	my ($tfname, $basedir, $app_tmpname) = @_;

	if(! -d $basedir) {
		#need an existing base (temporary) directory
		die "Directory $basedir does not exist.";
	}

	my $app_tmpdir = "$basedir/$app_tmpname/";
	if(! -d $app_tmpdir) {
		system("mkdir $app_tmpdir");
	}

	my $template = "$tfname.XXXXXXXXXX";
	my $fasta_tmpdir = tempdir ( $template, DIR => $app_tmpdir);
	return $fasta_tmpdir;
}



###########################################################################
# FASTA
###########################################################################
sub open_genome_file_keep_N
{
	my ($fn) = @_;

	open(IN, "< $fn") or die "Cannot open $fn for read\n";

	my @new_in = ();
	foreach my $s (<IN>) {
		$s =~ s/[\r\n]//g;
		if($s !~ /^\>/) { #if not tag
			push @new_in, $s;
		}
	}
	close(IN);

	my $in_join = join("", @new_in);

	my $name = File::Basename::fileparse($fn);
	my $seqlen = length($in_join);

	my %hash = (
		seq => $in_join,
		name => $name,
		filename => $fn,
		seqlen => $seqlen,
	);

	printf STDERR ("$fn genome file has %d characters\n", $seqlen);

	return \%hash;
}


sub open_genome_file
{
	my ($fn) = @_;

	open(IN, "< $fn") or die "Cannot open $fn for read\n";

	my @new_in = ();
	foreach my $s (<IN>) {
		chomp $s;
		$s =~ s/[\r\n]//g;
		if($s !~ /^\>/) { #if not tag
			push @new_in, $s;
		}
	}
	close(IN);

	my $in_join = join("", @new_in);
	$in_join =~ s/[^ACGTacgt]//g; #filter out all invalid characters
	my $seqlen = length($in_join);

	my $name = File::Basename::fileparse($fn);

	my %hash = (
		seq => $in_join,
		name => $name,
		filename => $fn,
		seqlen => $seqlen,
	);

	return \%hash;
}

sub open_genome_file_filter_coord {
	my ($fn, $coord_fn) = @_;

	my $genome = open_genome_file($fn);
	my $coord_lst = open_coord_file($coord_fn);
	my $seq = $genome->{seq};

	for(my $i = 0; $i < scalar(@$coord_lst); $i++) {
		my $offset = $coord_lst->[$i]->{offset};
		my $seqlen = $coord_lst->[$i]->{len};
		substr($seq, $offset, $seqlen) = 'X' x $seqlen;
	}
	my $len_old = length($seq);
	$seq =~ s/X//g;
	my $len_new = length($seq);
	my $len_diff = $len_old - $len_new;

	$genome->{seq} = $seq;

	print STDERR "Genome file filtered with coordinates: len_diff=$len_diff, len_old=$len_old, len_new=$len_new\n";
	print STDERR "  genome-file: $fn\n";
	print STDERR "  coord-file: $coord_fn\n";
	
	return $genome;
}

sub open_coord_file{
	my $file = shift;

	open(IN, "< $file") or die "Cannot open $file for read\n";

	my @lst = ();
	foreach my $ln (<IN>) {
		$ln =~ s/[\r\n]//g;
		$ln =~ s/\s//g;
		my @cols = split(/\,/,$ln);
		my %hash = (
			offset => $cols[0],
			len => $cols[1],
		);
		push(@lst, \%hash);
	}
	close(IN);
	return \@lst;
}


sub validate_coord
{
	my ($genome_fn, $coord_fn, $fsa_fn) = @_;
	
	my $genome_struct = open_genome_file($genome_fn);
	my $genome = $genome_struct->{seq};
	my $coord_lst = open_coord_file($coord_fn);
	my $fsa_lst = open_fasta_file($fsa_fn);
	if(scalar(@$fsa_lst) != scalar(@$coord_lst)) {
		die "Error: the number of seqs in FASTA and *.coord file are not the same";
	}

	for(my $i = 0; $i < scalar(@$fsa_lst); $i++) {
		my $offset = $coord_lst->[$i]->{offset};
		my $seqlen = $coord_lst->[$i]->{len};
		my @seq = split(//, $fsa_lst->[$i]);
		my @gsegment = split(//, uc(substr($genome, $offset, $seqlen)) );
		if(scalar(@seq) != $seqlen) {
			print STDERR "ERROR: validate_coord() found inconsistent seqlen for seq $i: $seqlen, " . scalar(@seq);
		}
		my $diff_count = 0;
		for(my $j = 0; $j < $seqlen; $j++) {
			if($seq[$j] ne $gsegment[$j]) {
				$diff_count++;
			}
		}
		if($diff_count > $seqlen / 4) {
			print STDERR "WARNING: coordinates inconsistent at validate_coord().\n";
			print STDERR "genome_fn: $genome_fn\n coord_fn: $coord_fn\n fsa_fn:$fsa_fn\n";
			print STDERR "seqind: $i, diff_count: $diff_count, seqlen: $seqlen.\n";
		}
	}
	print STDERR "Coordinates validation performed::\n";
	print STDERR " genome_fn: $genome_fn\n coord_fn: $coord_fn\n fsa_fn:$fsa_fn\n";
}

sub open_fasta_file {
	my $file = shift;

	open(IN, "< $file") or die "Cannot open $file for read\n";

	my @seqs = ();
	my $s = '';
	foreach my $ln (<IN>) {
		$ln =~ s/[\r\n]//g;
		$ln =~ s/\s//g;
		if($ln =~ /^\>/) { #if not tag
			if($s ne '') {
				push @seqs, $s;
				$s = '';
			}
		}
		else {
			$s .= uc($ln);
		}
	}
	if($s ne '') {
		push @seqs, $s; #add the last seq
	}
	close(IN);

	return \@seqs;
}

sub parseFa
{
	my $fname = shift;
	open(FILE, "< $fname") or die "Couldn't open file $fname: $!\n";
	my @istream = <FILE>;
	close(FILE);

	my $title = '';
	my $seq = '';
	my @lst = ();
	foreach my $s (@istream) {
		$s =~ s/[\r\n]//g;

		if(length($s) == 0) {
			next;
		}
		elsif($s =~ /^\>(.+)/) {
			my $newTitle = $1;
			if($seq ne '') {
				push(@lst, parseFa_createStruct($seq, $title));
				$seq = '';
			}
			$title = $newTitle;
		}
		else {
			$s =~ s/\s//g;
			$s =~ s/\*//g; #pir format has '*' that are useless
			$seq = $seq . $s;
		}

	}
	#if($seq =~ /[ACGTNXacgtnx\-\.]/) {
	if($seq =~ /\S+/) {
		push(@lst, parseFa_createStruct($seq, $title));
	}
	printf STDERR ("$fname has %d sequences.\n", scalar(@lst));
	return \@lst;
}

sub keepACGTonly
{
	my ($fsa) = @_;
	foreach my $s (@$fsa) {
		$s->{seq} =~ s/[^ACGTacgt]//g;
	}
	return $fsa;
}

sub isACGT
{
	my ($c) = @_;
	return ($c =~ /^[ACGTacgt]$/);
}

sub parseFa_createStruct
{
	my $seq = shift;
	my $title = shift;

	my %hash = (
		tag => $title,
		seq => uc($seq),
	);
	return \%hash;
}

sub writeFa
{
	my $fh = shift;
	my $fasta = shift;

	foreach my $f (@$fasta) {
		printf $fh (">%s\n", $f->{tag});
		printf $fh ("%s\n", $f->{seq});
	}
}

1;

__END__
