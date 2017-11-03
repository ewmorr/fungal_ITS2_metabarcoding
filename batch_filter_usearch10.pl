#/usr/bin/perl
#Eric Morrison
#2-9-15

use strict;
use warnings;

sub usage(){
	print STDERR q(

Usage: perl batch_filter_usearch.pl [input directory]

This script runs usearch10 expected error quality filtering on a directory of fastq files.

);
exit;
}
if( scalar(@ARGV) == 0 || $ARGV[0] eq "-h")
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

system "mkdir $workDir"."_usearchFilter";



foreach my $dir(@dirs)
	{
	system "mkdir $workDir"."_usearchFilter/$dir";
	system "usearch10 -fastq_filter $workDir/$dir/$dir.join.fastq -fastaout $workDir"."_usearchFilter/$dir/$dir.filter.fasta -fastqout $workDir"."_usearchFilter/$dir/$dir.filter.fastq -fastq_maxee 0.5 -relabel $dir. -threads 1";
	system 'sed "-es/^>\(.*\)/>\1;barcodelabel='."$dir".';/" < '.$workDir."_usearchFilter/$dir/$dir.filter.fasta > $workDir"."_usearchFilter/$dir/$dir.filter.label.fasta";
	#print "usearch10 -fastq_filter $workDir/$dir/$dir.join.fastq -fastaout $workDir"."_usearchFilter/$dir/$dir.filter.fasta -fastq_maxee 0.5 -relabel $dir. -threads 1\n";
	#print 'sed "-es/^>\(.*\)/>\1;barcodelabel='."$dir".';/" < '.$workDir."_usearchFilter/$dir/$dir.filter.fasta > $workDir"."_usearchFilter/$dir/$dir.filter.label.fasta\n";
	}

