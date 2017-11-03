#!/usr/bin/perl

use strict;
use warnings;

sub usage(){
	print STDERR q(
Usage: perl rmlt150.pl [input.fasta] [output_one.fasta] [output_two.fasta]

This script takes a fasta file as input and outputs sequences of greater than or equal to 150bp length to the first output file, and sequences less than 150bp lewngth to the second output file.

);
exit;
}
if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $input = $ARGV[0];
my $out = $ARGV[1];
my $out2 = $ARGV[2];

open (IN, "$input") || die "Can't open input\n";
open (OUT, ">$out") || die "Can't open output\n";
open (OUT2, ">$out2") || die "Can't open output\n";

chomp(my @fasta = <IN>);

my $fasta = join(":::::::::", @fasta);
my @seqs = split(">", $fasta);
shift @seqs;

my %fasta;
foreach my $seq(@seqs)
	{
	my @seq = split(":::::::::", $seq);
	my $head = shift @seq;
	$fasta{$head} = join("", @seq);
	}

my $tot = keys %fasta;
my $good;

foreach my $head(keys %fasta)
	{
#	print length($fasta{$head}), "\n";
	if (length($fasta{$head}) > 149)
		{
		$good++;
		print OUT ">",$head, "\n", $fasta{$head}, "\n";
		}
	elsif(length($fasta{$head}) < 150)
		{
		print OUT2 ">",$head, "\n", $fasta{$head}, "\n";
		}
	}
print "Total seqs: $tot\nGood seqs: $good\n";

