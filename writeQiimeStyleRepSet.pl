#!/usr/bin/perl
#Eric Morrison
#3-10-15
#This script rewrites rep_set and otus files from usearch to be compatible with qiime
use strict;
use warnings;

sub usage(){
	print STDERR q(
Usage: perl writeQiimeStyleRepSet.pl [rep_set.fasta] [otus.txt] [rep_set.qiimeStyle.fasta]

This script rewrites a usearch rep_set.fasta file to QIIME formatted fasta headers. The otu file is used to confirm that OTU numbers are accurate.

);
exit;
}
if( scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $rep = $ARGV[0];
my $otus = $ARGV[1];
my $repOut = $ARGV[2];

open(REP, "$rep") || die "Can't open rep set.\n";
open(OTUS, "$otus") || die "Can't open uc.\n";
open (OUT, ">$repOut") || die "Can't open rep_set output.\n";

chomp(my @rep = <REP>);
chomp(my @otus = <OTUS>);


my $repJ = join(":::::", @rep);
my @repS = split(">", $repJ);
shift(@repS);
my %rep;
my %size;#make hash of otu size indexed by sequence header to sort sequence hash by otu size. This should facilictate faster searching as the largest otus should be first in the otu list.
foreach my $repS(@repS)
	{
	my @repSS = split(":::::", $repS);
	my $head = shift(@repSS);
	$rep{$head} = join("", @repSS);
	$head =~ /.+;barcodelabel=.+;size=(\d+);/;
	$size{$head} = $1;
	}

foreach my $head(sort{ $size{$b} <=> $size{$a} } keys %size)
	{
	$head =~ /(.+\.\d+;barcodelabel=.+;)size=\d+;/;
	my $label = $1;
#	my $tog = 0;
	
	for(my $i = 0; $i < @otus; $i++)
		{
		my @otu = split("\t", $otus[$i]);
		my $otuNum = shift(@otu);
		
			if($label eq $otu[0])
				{
			if($otuNum == 642)
				{
				print "match\n";
				}
				$label =~ /(.+);barcodelabel=.+;/;
				my $newLabel = $1;
				$newLabel =~ s/\./_/;
				print OUT ">", $otuNum, " ", $newLabel, "\n", $rep{$head}, "\n";
#				$tog = 1;
				last;
				}
			
#		if($tog == 1)
#			{
#			splice(@otus, $i, 1);
#			last;
#			}
		}
	}

