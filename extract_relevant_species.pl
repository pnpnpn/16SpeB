#!/usr/bin/perl -w
use strict;
use FastaIO;

my $DEBUG0 = 1;

my @species_lst= (
	'Bacteroides fragilis',
	'Clostridium bifermentans',
	'Chlamydia trachomati',
	'Corynebacterium diphtheriae',
	'Haemophilus influenzae',
	'Helicobacter pylori',
	'Listeria monocytogenes',
	'Mycobacterium leprae',
	'Mycobacterium tuberculosis',
	'Mycoplasma hominis',
	'Neisseria gonorrhoeae',
	'Neisseria meningitidis',
	'Staphylococcus aureus',
	'Streptococcus pneumoniae',
	'Yersinia pestis',
);

main();



sub main {
	if(scalar(@ARGV) <= 1) {
		print "usage: <program> <destdir> [LIST OF FILES]\n";
		print "\n";
		print "./extract_relevant_species.pl relevant_species/ input_data/*.fasta\n";
		print "\n";
		exit(1);
	}

	print "$0 " . join(" ", @ARGV) . "\n\n";
	my $destdir = $ARGV[0];
	my @filelst = @ARGV[1 .. scalar(@ARGV)-1];

	#my @fsalst = ();
	#foreach my $f (@filelst) {
	#	push(@fsalst, parseFa($f));
	#}

	if(! -d $destdir) {
		die "Destination directory $destdir does not exist";
	}

	my %fsahash = ();
	foreach my $species (@species_lst) {
		$fsahash{$species} = [];
	}

	foreach my $file (@filelst) {
		extract_fasta_with_strlst(\%fsahash, \@species_lst, $file);
		print "\n";
	}

	foreach my $species (@species_lst) {
		(my $fileroot= $species) =~ s/\s+/\_/g;

		printf ("Total %5d sequences for $species\n", scalar(@{$fsahash{$species}}));
		my $fsafname = sprintf("$destdir/%s.fsa", $fileroot);
		open(my $outfh, ">$fsafname") or die "cannot open $fsafname for write: $!";
		writeFa($outfh, $fsahash{$species});
		close($outfh);
	}
}

sub extract_fasta_with_strlst {
	my ($fsahash, $strlst, $filename) = @_;

	my $fsa = parseFa($filename);

	foreach my $str (@$strlst) {
		my @new_seqset = ();
		foreach my $s (@$fsa) {
			if($s->{tag} =~ /$str/) {
				push(@new_seqset, $s);
			}
		}
		printf ("%5d sequences of $str from $filename\n", scalar(@new_seqset));

		push(@{$fsahash->{$str}}, @new_seqset);
	}
}


sub extract_fasta_with_str {
	my ($str, $filelst) = @_;

	#if(scalar(@$filelst) != scalar(@$fsalst)) {
	#	die "inconsistent fsalst and filelst";
	#}

	my @fsa_all = ();
	for(my $i = 0; $i < scalar(@$filelst); $i++) {
		my $filename = $filelst->[$i];
		my $fsa = parseFa($filename);

		my @new_seqset = ();
		foreach my $s (@$fsa) {
			if($s->{tag} =~ /$str/) {
				push(@new_seqset, $s);
			}
		}
		printf ("$str has %5d sequences from $filename\n", scalar(@new_seqset));

		push(@fsa_all, @new_seqset);
	}

	printf ("$str has %5d sequences total\n", scalar(@fsa_all));

	return \@fsa_all;
}


