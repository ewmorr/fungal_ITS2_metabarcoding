#!/usr/bin/perl
#Eric Morrison
#3-10-15

use strict;
use warnings;

sub usage(){
	print STDERR q(

Usage: perl writeQiimeStyleOtuFile.pl [otus.txt] [otus.qiimeStyle.txt]

This script rewrites sequence ID's from usearch format to QIIME format.

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

chomp(my @in = <IN>);

foreach my $otu(@in)
	{
	my @otu = split("\t", $otu);
	my $otuNum = shift@otu;
	print OUT $otuNum, "\t";
	foreach my $entry (@otu)
		{
		$entry =~ /(.+\.\d+);barcodelabel=.+;/;
		my $new = $1;
		$new =~ s/\./_/;
		print OUT $new, "\t";
		}
	print OUT "\n";
	}

