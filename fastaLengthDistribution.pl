#!/usr/bin/perl
#Eric Morrison
#10/8/14
#Usage: fastaLengthDistribution.pl [input] [output]

use strict;
use warnings;

sub usage(){
	print STDERR q(

Usage: perl fastaLengthDistribution.pl [input.fasta] [output.txt]

This script takes a fasta file as input and outputs a text file with the length distribution of sequences (bp v. number sequences).

);
exit;
}
if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $in = $ARGV[0];
my $out = $ARGV[1];

open(IN, "$in") || die "Can't open input.\n";
open(OUT, ">$out") || die "Can't open output.\n";

chomp(my @seqs = <IN>);
my $seqs = join(":::", @seqs);
my @splitSeqs = split(">", $seqs);
print $splitSeqs[0], "\n";
shift @splitSeqs;

my %seqs;
foreach my $seq(@splitSeqs)
	{
	my @seq = split(":::", $seq);
	my $head = shift(@seq);
	$seqs{$head} = join("", @seq);
	}


my %lenHash;
foreach my $head (keys %seqs)
	{
	$lenHash{length($seqs{$head})}++;
	}
	
print OUT "Length\tCount\n";	
foreach my $len (sort keys %lenHash)
	{
	print OUT $len, "\t", $lenHash{$len}, "\n";
	}
