#!/usr/bin.perl
#Eric Morrison
#2-10-15
#This script creates a primer list for use with trimmomatic given a list associating
#SampleID with forward and reverse barcodes.

use strict;
use warnings;

sub usage(){
	print STDERR q
	(
Usage: perl trimmomaticPrimerBuild.pl [input primer list] [output directory]

This script takes a list of samples associated with barcodes and build primers lists for use with trimmomatic. Primers for each sample are output to a separate directory that is specified in the call to the script.

);
exit;
}
if((@ARGV==0) || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $in = $ARGV[0];
my $outDir = $ARGV[1];
system "mkdir $outDir";

open(IN, "$in") || die "Can't open barcode file.\n";


chomp(my @bc = <IN>);

if(scalar(@bc) == 1)
	{
	$bc[0] =~ s/\r|\n|\r\n/\n/g;
	@bc = split("\n", $bc[0]);
	}
	
foreach my $sample(@bc)
	{
	if($sample =~ /#.*/)
		{
		next;
		}
	my @sample = split("\t", $sample);
	open(TEMP, ">$outDir/$sample[0].fa") || die "Can't open ouput.\n";
	print TEMP ">Prefix".$sample[0]."/1\n",#We write two sets of pairs to account for degneracy in forward primer which is indicated by "A" or "G"
		"AATGATACGGCGACCACCGAGATCTACAC".$sample[1]."ACACTCTTTCCCTACACGACGCTCTTCCGATCT\n",#forward primer one
		">Prefix".$sample[0]."/2\n",
		"CAAGCAGAAGACGGCATACGAGAT".$sample[2]."GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n";#reverse primer one
		#">Prefix".$sample[0]."G/1\n",
		#"AATGATACGGCGACCACCGAGATCTACAC".$sample[1]."GCAGCGAGCC"."GG"."GTGAGTCATCGAATCTTTG\n",#forward primer two
		#">Prefix".$sample[0]."G/2\n",
		#"CAAGCAGAAGACGGCATACGAGAT".$sample[2]."GGTCTGCGCG"."AA"."TCCTCCGCTTATTGATATGC\n";#reverse primer two
	close(TEMP);
	}
