#/usr/bin/perl
#Eric Morrison
#2-9-15

use strict;
use warnings;

sub usage(){
	print STDERR q(
Usage: perl batch_count_fastq_entries.pl [input directory] [output file]
This script counts number of entries in the fastq files in a directory.

);
exit;
}

if(scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
	{
	&usage();
	}

my $workDir = $ARGV[0];
$workDir =~ s/\///;

my $out = $ARGV[1];
open(OUT, ">$out") || die "Can't open output.\n";

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


foreach my $dir(keys %files)
	{
	print OUT $dir, "\t";
	
	foreach my $name (@{ $files{$dir} })
		{
		system "expr \$(cat $workDir/$dir/$name | wc -l) / 4 > tempCount.txt";
		open(COUNT, "tempCount.txt") || die "Can't open count.\n";
#		print "expr \$(cat $workDir/$dir/$name | wc -l) / 4 > tempCount.txt\n";
		chomp(my $count = <COUNT>);
		system "rm tempCount.txt";
		close(COUNT);
		print OUT $name, "\t", $count, "\t";
		}
	print OUT "\n";
	}

