#!/usr/bin.perl
#Eric Morrison
#2-11-15
#This runs trimmomatic on paired Illumina reads when given a directory of reads and a directory of adapters.

use strict;
use warnings;

sub usage(){
	print q(
Usage: perl batchTrimmomatic.pl [input directory] [adaptors directory]

This script runs trimmomatic on files of paired Illumina reads and returns output to the specified directory.

);
exit;
}
if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $fastqDir = $ARGV[0];
my $adapDir = $ARGV[1];

$fastqDir =~ s/\///g;
$adapDir =~ s/\///g;

system "ls $fastqDir > tempFileList.txt";
my $files = "tempFileList.txt";
open(FILES, "$files") || die "Cant open file list.\n";
system "rm tempFileList.txt";
chomp(my @files =<FILES>);

system "ls $adapDir > tempPrimerList.txt";
my $primers = "tempPrimerList.txt";
open(PRIMER, "$primers") || die "Cant open primer list.\n";
#system "rm tempPrimerList.txt";
chomp(my @primers =<PRIMER>);

my %filesList;
my %primerList;
foreach my $name (@files)
	{
	$name =~ /(.*)(R1|R2)(.*)/;
	$filesList{$1} = $3;
	
	}
	
foreach my $name(keys %filesList)
	{
	$name =~ /^(.+)?_S\d+.*/;
	my $tempName = $1;
#	print $tempName, "\n";
	
	for(my $i = 0; $i<@primers; $i++)
		{
		my $tempPrim = "dummy";
		if(defined($primers[$i]) == 1)
			{
			$primers[$i] =~ /(.*)\.fa/;
			$tempPrim = $1;
			if($tempName eq $tempPrim)
				{
				$primerList{$name} = $tempPrim;
#				print $primerList{$name}, "\n";
				}
			}
		
		}
	}

system "mkdir $fastqDir"."_trim";
foreach my $name(keys %filesList)
	{
	system "mkdir $fastqDir"."_trim/$primerList{$name}";
	system "java -jar bin/trimmomatic-0.33.jar PE -phred33 -trimlog $fastqDir"."_trim/$primerList{$name}/trimLogFile.txt $fastqDir/$name"."R1"."$filesList{$name} $fastqDir/$name"."R2"."$filesList{$name} $fastqDir"."_trim/$primerList{$name}/$name"."R1"."pair.fastq $fastqDir"."_trim/$primerList{$name}/$name"."R1"."unpair.fastq $fastqDir"."_trim/$primerList{$name}/$name"."R2"."pair.fastq $fastqDir"."_trim/$primerList{$name}/$name"."R2"."unpair.fastq ILLUMINACLIP:$adapDir/$primerList{$name}.fa:1:30:15:1:true";
	#Write stout to file
	print "java -jar bin/trimmomatic-0.33.jar PE -phred33 -trimlog $fastqDir"."trim/$primerList{$name}/trimLogFile.txt $fastqDir/$name"."R1"."$filesList{$name} $fastqDir/$name"."R2"."$filesList{$name} $fastqDir"."trim/$primerList{$name}/$name"."R1"."pair.fastq $fastqDir"."trim/$primerList{$name}/$name"."R1"."unpair.fastq $fastqDir"."trim/$primerList{$name}/$name"."R2"."pair.fastq $fastqDir"."trim/$primerList{$name}/$name"."R2"."unpair.fastq ILLUMINACLIP:$adapDir/$primerList{$name}.fa:1:30:15:1:true\n";
#	open(TEMP, ">$outDir/$sample[0].fa") || die "Can't open ouput.\n";
	
	}
