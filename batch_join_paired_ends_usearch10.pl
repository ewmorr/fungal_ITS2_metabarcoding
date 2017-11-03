#/usr/bin/perl
#Eric Morrison
#2-9-15

use strict;
use warnings;

sub usage(){
	print STDERR q(
Usage: batch_join_paired_usearch10.pl [input directory]

This script runs userach10 read merging on a directory of files.

);
exit;
}

if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $workDir = $ARGV[0];
$workDir =~ s/\///;

system "ls $workDir > tempDirList.txt";
my $dirs = "tempDirList.txt";
open(DIRS, "$dirs") || die "Cant open directory list.\n";
system "rm tempDirList.txt";
chomp(my @dirs =<DIRS>);

my %files;
foreach my $dir(@dirs)
	{
	system "ls $workDir/$dir > tempFileList.txt";
	my $files = "tempFileList.txt";
	open(FILES, "$files") || die "Cant open file list.\n";
	system "rm tempFileList.txt";
	chomp(my @files =<FILES>);
	$files{$dir} = [ @files ];
	}

my %filesList;
foreach my $dir(keys %files)
	{
	foreach my $name (@{ $files{$dir} })
		{
		$name =~ /(.*)(R1|R2)pair.fastq/;
		$filesList{$dir} = $1;
		}
	}

system "mkdir $workDir"."_merged";
foreach my $name(keys %filesList)
	{
	system "mkdir $workDir"."_merged/$name";
	system "usearch10 -fastq_mergepairs $workDir/$name/$filesList{$name}"."R1"."pair.fastq -reverse $workDir/$name/$filesList{$name}"."R2"."pair.fastq -fastqout $workDir"."_merged/$name/$name.join.fastq -fastqout_notmerged_fwd $workDir"."_merged/$name/$name.un1.fastq -fastqout_notmerged_rev $workDir"."_merged/$name/$name.un2.fastq -fastq_minovlen 20 -minhsp 20 -fastq_minqual 3 -fastq_maxdiffs 8";
	}
