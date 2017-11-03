#!/usr/bin/perl
# getAllSeqsMothurOTUs.pl
# Eric Morrison
# 6/08/2017
# Usage: getAllSeqsMothurOTUs.pl [cluster_results.uc] [derep.uc] [ouput.txt]
#
# This script reads a uc formatted file from usearch and
# and a derep uc file.
# Returns a Mothur formatted file with headers for each 
# sequence in an OTU. 


use strict;
use warnings;

sub usage(){
	print STDERR q(
Usage: perl getAllSeqsMothurOTUs.pl [cluster_results.uc] [derep.uc] [output.txt]

This script reads two usearch .uc files, the first from a final clustering, and the second from dereplication. The script searches sequence id's in each file to find the sequences that belong to OTUs and print OTUs in a mothur format OTU file.

);
exit;
}
if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $clust = $ARGV[0];
my $derep = $ARGV[1];
my $out = $ARGV[2];

open (CLUST, "$clust") || die "Can't open input file.\n";
open (DEREP, "$derep") || die "Can't open fasta file.\n";
open (OUT, ">$out") || die "Can't open output file.\n";

chomp(my @clust = <CLUST>);
chomp(my @derep = <DEREP>);

my %clusters;
my %derep;
my %otus;

# Hash for derep seqs. Each cluster is placed in an array indexed in hash by seed sequence header.
#Seeds are not included in array as these are in primary cluster data
my $count1;
my $count2;
foreach my $line (@derep)
	{
	my @derepLine = split("\t", $line);
	
	if ($derepLine[0] eq "H")
		{
		$count1++;
		push @{ $derep{$derepLine[9]} } , $derepLine[8];
		}
#Do not need to store seeds of derep clusters. These will be sequences to search primary cluster array
#but are already included in array by definition
	if ($derepLine[0] eq "S")
		{
		$count2++;
#		unshift @{ $derep{$derepLine[8]} }, $derepLine[8];
		}
	}
my $count = $count1 + $count2;
print $count, " sequences being analyzed. $count2 unique\n";


# Hash for sequence clusters. Each cluster is place in an array indexed in hash by 
# cluster number. Must split header by ; to remove sequence count data. Seed is array
# index [0] to facilitate search by derep seeds as cluster seeds are more likely to be
# large sequence clusters.
foreach my $line (@clust)
	{
	my @clustLine = split("\t", $line);
	if ($clustLine[0] eq "H")
		{
		my @splitHead = split(";", $clustLine[8]);
		push @{ $clusters{$clustLine[1]} } , $splitHead[0].";".$splitHead[1].";";
		}
	if ($clustLine[0] eq "S")
		{	
		my @splitHead = split(";", $clustLine[8]);
		unshift @{ $clusters{$clustLine[1]} } , $splitHead[0].";".$splitHead[1].";";
			
		}
	}	

# Search derep hash for clusters that match those in sequence clusters.
# Write a new hash indexed by OTU number with all 
# sequence headers in each cluster. 
my %finalClust;		
foreach my $otuNum(keys  %clusters)
	{
	$finalClust{$otuNum} = [@{ $clusters{$otuNum} }];
#	print $otuNum, "\n";
	foreach my $clust (@{ $clusters{$otuNum} })	
		{
		my $searchNum = $clust;
	#	print $searchNum, "\n";
		if(defined($derep{$searchNum}) == 1) 
			{
			push (@{ $finalClust{$otuNum} }, @{ $derep{$searchNum} } );
			}
		}
	}
# Print clusters to file

foreach my $final(sort {$a <=> $b} keys %finalClust)
	{
	my $otuID = $final;
	$final =~ s/q//;
	print OUT "$final\t";
	foreach my $id (@{ $finalClust{$otuID} })
		{
		print OUT $id, "\t";
		}
	print OUT "\n";
	}

