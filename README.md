A bioinformatic pipeline for processing fungal ITS sequence data from MiSeq dual-index metabarcoding sequencing projects

Author: Eric W. Morrison, 7/1/2015, revised 10/2/2015

Description: Pipeline for processing MiSeq paired-end read data for fungal meta-barcoding experiments. Data is assumed to be demultiplexed by the sequencing facility. Batch processing scripts are provided to run processes on directories of demultiplexed sequence files.


##########################################################################################
										###NOTE###

This document is based on a tutorial I put together for novices to metabarcoding. The tutorial used a set of example files to work through the basic process using data from one sample (example_sample_1), or a set of two files to be processed together. 
I have not made much attempt to remove the "tutorial-esque" language in this document, so please excuse if it reads at a basic level.

##########################################################################################



Required software:\n
For removal of MiSeq adaptors\n
Trimmomatic-0.33 <http://www.usadellab.org/cms/?page=trimmomatic>

For read merging, quality filtering, sequence dereplicating, and clustering
Usearch v. 10.0 <http://www.drive5.com/usearch/>
###or latest version. usearch is updated frequently and this code was last tested with version 10, and originally written with v. 8 in mind. Some commands may need update.

For ITS region extraction
ITSx  v. 1.0.7 <http://microbiology.se/software/itsx/>

Perl scripts
trimmomaticPrimerBuild.pl
batchTrimmomatic.pl
batch_count_fastq_entries.pl
fastaLengthDistribtuion.pl
batch_join_paired_ends_usearch8.pl
batch_filter_usearch.pl
rmlt150.pl
combine_derep_and_clustering_uc_files.pl
writeQiimeStyleRepSet.pl
writeQiimeStyleOtuFile.pl

Additionally:
We’ll need the UNITE dynamically clustered reference database and a working install of QIIME if using the blast wrapper and OTU table building function there is desired. 

Pipeline flow:
Step 1. Remove adaptors (Trimmomatic, trimmomaticPrimerBuild.pl, batchTrimmomatic.pl)
Step 2. Merge paired end reads (usearch8, batch_join_paired_ends_usearch8.pl) 
Step 3. Quality filter reads (usearch8, batch_filter_usearch.pl)
Step 4. Filter by length; recommended minimum length 150bp (rmlt150.pl)
Step 5. Dereplicate sequences (usearch8)
Step 6. Extract ITS region (ITSx)
Step 7. Cluster sequences; recommended 97% sequence similarity (usearch8)
Step 8. Build OTU map file for input into QIIME or other software (getAllSeqsMothurOTUs.pl)
Step 9. Assign taxonomy and make OTU table (QIIME)

Note: The original sequence files used for input to this pipeline are assumed to be in gzipped fastq or fastq format (file extension .fastq.gz or .fastq). For batch processing it is assumed that file names include the sample identifiers associated with each barcode, and include the designation “R1” or “R2” for reverse or forward reads (this is standard for paired-end MiSeq data). Also, IT IS ALWAYS A GOOD IDEA TO CHECK THE NUMBER OF SEQUENCES REMAINING AFTER PROCESSING STEPS, especially steps that involve quality filtering, to see that the number of sequences retained matches expectations.  For steps that may affect the length distribution of sequences (e.g. read merging) sequence length distribution should also be examined.


For the purpose of this tutorial we will work with the sequence data files for samples “example_sample_1” and “example_sample_2.” We start with the example_sample_1 read files in a directory called “miseq_test” and the sequence files for both samples in a directory called “miseq_test_batch.” The sample barcode list is provided in the file “miseq_test_primers.txt.” We will enter the commands to process the single sample data files for illustration purposes, but most experiments will have many samples that need to be processed. The batch processing scripts are designed to run each process on multiple sample files so that the user is saved tedious typing etc. 


Step 1. Remove adaptors
This step removes sequencing adaptors and PCR primers from sequence reads. Sample barcodes can interfere with down-stream processing (especially read merging, ITS extraction, and clustering) as barcodes and adaptors are arbitrary sequences with respect to the biological sequence information. Note that the sequences are demultiplexed coming off of the sequencer, but there may be cases where we sequenced through the primers at the opposite end of the read, or primers remain for other reasons.

In order to run Trimmomatic we need to provide the primer sequences we are searching for. Putting together a separate fasta file of primers for each sample can be quite time consuming (and error prone), so we use a perl script.

```
#Build primers
perl scripts/trimmomaticPrimerBuild.pl miseq_test_primers.txt miseq_test_primers_dir/

#Look at output
ls miseq_test_primers_dir/
more miseq_test_primers_dir/example_sample_1.fa
```

Note that many perl scripts are written to give usage information when the script is run without arguments like so:

```
perl scripts/trimmomaticPrimerBuild.pl 

        Usage: perl trimmomaticPrimerBuild.pl [input primer list] [output directory]

        This script takes a list of samples associated with barcodes and build primers lists for use with trimmomatic. Primers for each sample are output to a separate directory that is specified in the call to the script.
```

The script creates a directory of adaptor files. The input is a tab-delimited text file containing a list of barcodes associated with sample IDs in the following format:

```
#SampleID	Index2sequence	Index1sequence
samp1	AAGCAGCA	ACCTAGTA
samp2	ACGCGTGA	ACCTAGTA
samp3	CGATCTAC	ACCTAGTA
samp4	TGCGTCAC	ACCTAGTA
samp5	GTCTAGTG	ACCTAGTA
samp6	CTAGTATG	ACCTAGTA
samp7	GATAGCGT	ACCTAGTA
samp8	TCTACACT	ACCTAGTA
```

The sample ID’s listed should be the same as those indicated in the fastq file names for trimming. Note that Index1 is the barcode attached to the reverse primer (this barcode is sequenced before Index2 in the sequencing reaction). If only one barcode was used the appropriate column can be left blank and/or the script modified for input of one barcode. The pad, linker, and primer sequences are hardcoded in the script and these will need to be modified to the appropriate sequence for the experiment. Alternatively, one could build these files by hand according to the specifications given in the Trimmomatic manual.

In order to run adaptor trimming with Trimmomatic for a single set of paired-end read files run the following command in the command line:

```
#run trimmomatic on a pair of files

mkdir miseq_test_trim

mkdir miseq_test_trim/example_sample_1

java -jar bin/trimmomatic-0.33.jar PE -phred33 -trimlog miseq_test_trim/example_sample_1/trimLogFile.txt miseq_test/example_sample_1_S10_L001_R1_001.fastq.gz miseq_test/example_sample_1_S10_L001_R2_001.fastq.gz miseq_test_trim/example_sample_1/example_sample_1_S10_L001_R1pair.fastq miseq_test_trim/example_sample_1/example_sample_1_S10_L001_R1unpair.fastq miseq_test_trim/example_sample_1/example_sample_1_S10_L001_R2pair.fastq miseq_test_trim/example_sample_1/example_sample_1_S10_L001_R2unpair.fastq ILLUMINACLIP:miseq_test_primers_dir/example_sample_1.fa:1:30:15:1:true
```

In most cases all of the reads should end up in the paired read one and paired read two sequence files. Any sequences in the unpaired files are those that did not have a mate, and may indicate a problem with the sequence or with the data file. 

The options specified after the ILUMMINACLIP argument determine how the program will proceed with adaptor trimming. The Trimmomatic manual gives a detailed explanation of these settings. The adaptors file format is essential for proper adaptor trimming. See the Trimmomatic manual for a description of the file format.

Adaptor trimming can be run on a directory of fastq files using the script batchTrimmomatic.pl, which runs Trimmomatic palindromic adaptor trimming on each set of paired read one/read two files. A list of primers associated with each set of sample files is required. One separate primer file for each sample pair should be provided in a single directory that is separate from the directory containing fastq files. Output will be written to a new directory that is named the same as the input directory but with “_trim” appended (e.g. input_fastq_directory_trim). One subdirectory will be generated within this directory for each sample with the sample ID as directory name, and the output files for each sample written to the sample subdirectory. 

```
#run trimmomatic on a directory of files
perl scripts/batchTrimmomatic.pl miseq_test_batch/ miseq_test_primers_dir/
```

Checking numbers of sequences in fastq files:
The number of sequences in each fastq file can be counted using ‘cat’:

```
expr $(cat miseq_test_batch_trim/example_sample_1/example_sample_1_S10_L001_R1pair.fastq | wc -l) / 4
expr $(cat miseq_test_batch_trim/example_sample_1/example_sample_1_S10_L001_R1unpair.fastq | wc -l) / 4
```

This counts the total number of lines in the file and divides by four. There are four lines per sequence in standard fastq format, though there are exceptions to this, and unexpected results at this stage may be due to inconsistent file formats.

The number of sequences in each file in a directory of fastq files can be checked using the script batch_count_fastq_entries.pl:

```
#Count fastq entries in a directory of files
perl scripts/batch_count_fastq_entries.pl miseq_test_batch_trim/  trimmed_counts.txt
less trimmed_counts.txt
```

Step 2. Merge reads

There are several programs available for read merging; we will use usearch (Edgar 2013) here. The usearch author provides some strong rationale for his approach (Edgar & Flyvbjerg 2015, http://drive5.com/usearch/manual/merge_pair.html, http://drive5.com/usearch/manual/quality_score.html), and we use the program in other parts of the pipeline. The command for a pair of files:

```
#Merge reads
#Run usearch10 read merging on a pair of files
mkdir miseq_test_trim_merged
mkdir miseq_test_trim_merged/example_sample_1
usearch10 -fastq_mergepairs miseq_test_trim/example_sample_1/example_sample_1_S10_L001_R1pair.fastq -reverse miseq_test_trim/example_sample_1/example_sample_1_S10_L001_R2pair.fastq -fastqout miseq_test_trim_merged/example_sample_1/example_sample_1.join.fastq -fastqout_notmerged_fwd miseq_test_trim_merged/example_sample_1/example_sample_1.un1.fastq -fastqout_notmerged_rev miseq_test_trim_merged/example_sample_1/example_sample_1.un2.fastq -fastq_minovlen 20 -minhsp 20 -fastq_minqual 3 -fastq_allowmergestagger –fastq_maxdiffs 8
```

A description of these parameters is given in the usearch manual (https://www.drive5.com/usearch/manual/cmd_fastq_mergepairs.html). The critical parameters in this case are the ‘-fastq_minovlen’ argument, which specifies the minimum overlap required to merge a pair of reads, and ‘-fastq_maxdiffs,’ which specifies the total allowable mismatches. We set –fastq_minovlen value to 20 as we expect our read length distribution for ITS2 reads to have a maximum of 450 to 500bp. Assuming 250bp reads, a minimum overlap of 16 leaves us a maximum read length of 250 * 2 - 20 = 480 for merged reads, so that we should not lose too many reads due to technical bias (although we may lose some). Note that the number of possible combinations of nucleotides for a sequence of length n is 4^n, so for an overlap of 20 nucleotides and allowing for eight mismatches (i.e. wildcard nucleotides in the overlap) we have 4^12 = 16,777,216 possible sequences. The minimum overlap can be adjusted if needed, but note that the ‘-minhsp’ must be at least as low as the minimum overlap. The usearch approach attempts to correct mismatched bases in the alignment according to their quality scores (Edgar & Flyvbjerg 2015). So, setting –fastq_maxdiffs at too low a value may cause us to discard useful reads. Setting too high a value (e.g. 100%), could potentially result in spurious read alignment. At this stage some examination of the affects of different parameters and comparison of those results to reasonable expectations is probably helpful. Note that we also set the minimum per base quality score to 2 using ‘-fastq_minqual 3’ as very low quality base calls can interfere with read merging.

To run the batch processing script using the parameters described above:

```
#on a directory
perl scripts/batch_join_paired_ends_usearch8.pl miseq_test_batch_trim/
```

The output will be written to a new directory with “_merged” appended (e.g. input_directory_merged). Note that in this case the directory for input is a directory containing one subdirectory for each sample (as is generated by batchTrimmomatic.pl).
The number of sequences retained should be checked after this step and compared to the number of input sequences (i.e. after adaptor trimming). A large loss of sequences at this stage may indicate a problem with the merge parameters or input files, or low sequence quality. It is also a good idea to check the read length distribution at this stage to see that this meets expectations. For example, Fig. 1 shows read length distribution of fungal ITS2 sequences from Illumina MiSeq (a) and 454 sequencing (b) from similar environments. The 454 data does not require read merging and we would expect read lengths from the different platforms to be roughly similar at this stage.

Read length distribution can be checked using the -fastq_stats command in usearch10:

```
#Sequence read lengths single file
usearch10 -fastq_eestats2 miseq_test_trim_merged/example_sample_1/example_sample_1.join.fastq
```

For a directory of files we can concatenate the separate files to a single file:

```
#Sequence read lengths directory of files
cat miseq_test_batch_trim_merged/*/*.join.fastq > merged_batch.fastq
usearch10 -fastq_eestats2 merged_batch.fastq
```

Step 3. Quality filtering
Some rationale for our approach:
There are several different approaches/philosophies concerning quality filtering of metabarcoding data:

Bokulich et al. (2013) Quality filtering vastly improves diversity estimates from Illumina amplicon sequencing. Nature Methods.
Kozich et al. (2013) Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. Applied and Environmental Microbiology.
Edgar (2013) UPARSE: highly accurate OTU sequences from microbial amplicon reads. Nature Methods.
Edgar & Flyvbjerg (2015) Error filtering, pair assembly and error correction for next-generation sequencing. Bioinformatics.

However, most approaches call for relatively stringent quality screening. Stringent quality filtering is often not a concern in genome sequencing of pure cultures, for example, because high read depth allows filtering of erroneous reads by taking a consensus of reads for an individual base. That is, different genes from an individual look different enough in most cases that a few errors will not affect our interpretation of which reads belong where, and we have many reads for each region allowing us to determine the correct read by consensus (i.e. the most common sequence at a site), especially with help from long-read sequencing technology. This is not the case in metabarcoding studies where we sequence a single (usually small) gene or region from many organisms, and attempt to use sequence similarity to define members of a community. In this case we cannot use consensus to find the correct sequence as we expect many sequences to be relatively similar, and we must trust base calls to define our operational taxonomic units. In addition, we expect error rates from metabarcoding studies to be higher than those in genome sequence studies because of low sequence diversity, which reduces read quality in Illumina sequencing. Therefore, proper (strict) quality filtering is essential, and a lack of this may result in erroneously high richness estimates, type II error in assessment of treatment effects on community composition (assuming errors are randomly distributed among treatments), and high numbers of unidentifiable sequences (see Edgar 2013 for a discussion of the pitfalls of poor sequence quality screening).

Here we use the usearch “expected error” filtering approach with an expected error rate specification of 0.5 (i.e. probability of less than one (0.5) errors per read; http://www.drive5.com/usearch/manual/exp_errs.html, http://www.drive5.com/usearch/manual/expected_errors.html).

```
#Single file
mkdir miseq_test_trim_merged_usearchFilter
mkdir miseq_test_trim_merged_usearchFilter/example_sample_1
usearch10 -fastq_filter miseq_test_trim_merged/example_sample_1/example_sample_1.join.fastq -fastaout miseq_test_trim_merged_usearchFilter/example_sample_1/example_sample_1.filter.fasta -fastq_maxee 0.5 -relabel example_sample_1. -threads 1
```

For downstream steps we need to relabel the sequences in the newly produced fasta file:

```
sed "-es/^>\(.*\)/>\1;barcodelabel=example_sample_1;/" < miseq_test_trim_merged_usearchFilter/example_sample_1/example_sample_1.filter.fasta > miseq_test_trim_merged_usearchFilter/example_sample_1/example_sample_1.filter.label.fasta
```

Look at the output:

#Look at output
~$ head miseq_test_trim_merged_usearchFilter/example_sample_1/example_sample_1.filter.fasta
~$ head miseq_test_trim_merged_usearchFilter/example_sample_1/example_sample_1.filter.label.fasta

These are now in fasta format. To count the sequences:

```
#Count sequences
grep ">" miseq_test_trim_merged_usearchFilter/example_sample_1/example_sample_1.filter.label.fasta | wc
```

For length distribution we can use a perl script:

```
#Calculate length distribution
perl scripts/fastaLengthDistribution.pl miseq_test_trim_merged_usearchFilter/example_sample_1/example_sample_1.filter.label.fasta filtered_len_dist.txt
head filtered_len_dist.txt
tail filtered_len_dist.txt
```

For batch processing:

```
#For batch processing
~$ perl scripts/batch_filter_usearch10.pl miseq_test_batch_trim_merged/
```

The batch processing script writes results to a directory with _usearchFilter appended and automatically performs the relabeling step (results file .label.fasta). At this point we can concatenate the separate sample files to a single file. All downstream steps will proceed with the samples combined:

```
#Cat to a single file
cat miseq_test_batch_trim_merged_usearchFilter/*/*.label.fasta > miseq_test_batch_trim_merged_usearchFilter/all_seqs.fasta
```

Number sequences remaining and sequence length distribution can also be checked at this stage.

```
#Count sequences
grep ">" miseq_test_batch_trim_merged_usearchFilter/all_seqs.fasta | wc

#Length distribution
perl scripts/fastaLengthDistribution.pl miseq_test_batch_trim_merged_usearchFilter/all_seqs.fasta batch_filtered_len_dist.txt
head batch_filtered_len_dist.txt
tail batch_filtered_len_dist.txt
```

Step 4. Filter by sequence length
**From this step we will proceed with the concatenated sample file only, rather than processing both the single sample and the concatenated file.
 
It is common practice in metabarcoding studies to remove sequences below a certain length threshold based on the expected length of the biological sequence, and the expected length output of the sequencer. Here we use a minimum sequence length cutoff of 150 bp:

```
#Filter by sequence length
#Proceeding with only one file as the separate sample files are now concatenated

perl scripts/rmlt150.pl miseq_test_batch_trim_merged_usearchFilter/all_seqs.fasta miseq_test_batch_trim_merged_usearchFilter/seqs_ge_150.fasta miseq_test_batch_trim_merged_usearchFilter/seqs_lt_150.fasta
```

The first output file contains the sequences with length greater than or equal to 150bp and the second file contains sequences shorter than 150bp. The script prints the number of sequences analyzed (“total” sequences) and the number retained (“good” sequences), but the number of sequences in each file may also be checked with ‘grep’:

```
grep ">" miseq_test_batch_trim_merged_usearchFilter/seqs_ge_150.fasta | wc
grep ">" miseq_test_batch_trim_merged_usearchFilter/seqs_lt_150.fasta | wc
```

Step 5. Dereplicate sequences
At this stage we perform our first sequence-clustering step. The goal here is to combine all of the sequences that are exactly the same, resulting in a fasta file containing all of the unique sequences, as well as a map (“uc” file in usearch parlance) that describes which sequences were matched. One desirable outcome from dereplicating is that we significantly reduce the size of the dataset for downstream processing. This will be important for the next step, ITS extraction, which is demanding in terms of computer time. To dereplicate with usearch:

```
#Dereplicate
usearch10 -fastx_uniques miseq_test_batch_trim_merged_usearchFilter/seqs_ge_150.fasta -uc miseq_test_batch_trim_merged_usearchFilter/seqs_ge_150.derep.uc -fastaout miseq_test_batch_trim_merged_usearchFilter/seqs_ge_150.derep.fasta -sizeout
#Count unique sequences
grep ">" miseq_test_batch_trim_merged_usearchFilter/seqs_ge_150.derep.fasta | wc
```

The fasta file output will be used as input for the ITS extraction step, while the uc file will be used later to assemble an OTU map.

Step 6. ITS extraction
This step involves using Hidden Markov Model (HMM) searches to find sequences with homology to known 5.8S and LSU sequences, then trim these returning the ITS2 region. Input is the dereplicated sequence file and output is written to the specified directory.

```
#ITS extraction
#This is probably best to be run on a server if available. If run on a server use 'nohup' before the command to prevent the shell from sending a hangup signal as well as & after to run in background.
nohup ITSx -i miseq_test_batch_trim_merged_usearchFilter/seqs_ge_150.derep.fasta -o miseq_test_batch_ITSx --preserve T --reset T --cpu 8 -t F -E 0.001 &

#Enter to return to the prompt


#To monitor progress
top
less nohup.out


#When finished look at run summary
less miseq_test_batch_ITSx.summary.txt

#count seqs in output
grep ">" miseq_test_batch_ITSx.ITS2.fasta | wc
grep ">" miseq_test_batch_ITSx.ITS1.fasta | wc
grep ">" miseq_test_batch_ITSx_no_detections.fasta | wc

#move files to a directory
mkdir miseq_test_batch_ITSx/
mv miseq_test_batch_ITSx* ./miseq_test_batch_ITSx/
```

The ‘--cpu’ flag sets the number of CPU cores used by the program. On machines without multiple cores this should be left out of the command, but runtime will be increased significantly. The ‘-E’ flag sets the E-value cutoff for identifying a conserved sequence region (i.e. 5.8S and LSU in our case). The default setting is 10^-5 (lower values are more stringent), and running with 0.001 allows us to retain more sequences at the expense of potentially including false positive matches. Assuming we have a done a good job quality filtering and intend to remove unwanted taxa (e.g. plants) after taxonomic identification this should not be too much of a problem. Additionally, running with a stringent E-value cutoff may artificially reduce the dataset if we assume our dataset includes some previously unobserved sequences that may not have a good 5.8S or LSU match in the HMM search database. The user should make their own decision here based on the goals of the study, and philosophy regarding “unknowns.” It is important to check the number of sequences retained at this stage (in the .ITS2.fasta file) as large numbers of sequences can be discarded if parameters are set improperly.

Step 7. OTU clustering

We will now cluster the extracted ITS2 sequences into OTUs based on sequence similarity. This is done in usearch v. 8 using several steps. First we sort by cluster size to allow clustering in usearch:

```
#make a new directory for clustering output
~$ mkdir miseq_test_batch_cluster

#sort by cluster size
usearch10 -sortbysize miseq_test_batch_ITSx/miseq_test_batch_ITSx.ITS2.fasta -fastaout miseq_test_batch_cluster/seqs.sorted.fasta

#Look at ouput
head miseq_test_batch_cluster/seqs.sorted.fasta
```

Next we perform sequence similarity clustering:

```
#Find representative sequences and identify chimeras
usearch10 -cluster_otus miseq_test_batch_cluster/seqs.sorted.fasta -otus miseq_test_batch_cluster/rep_set.fasta -uparseout miseq_test_batch_cluster/cluster_results.txt -otu_radius_percent

#Count
grep ">" miseq_test_batch_cluster/rep_set.fasta | wc
```

The default sequence similarity threshold in usearch v. 10 is 97%. (Although see here <https://www.drive5.com/usearch/manual/cmd_cluster_otus.html> for arguments against 97% clustering). Note that there is no sequence map output (uc file) to tell us which sequences belong to which OTUs. This is a change from previous versions of usearch and requires running usearch’s “mapping” function:

```
#Map sequences to representative sequences
usearch10 -usearch_global miseq_test_batch_cluster/seqs.sorted.fasta -db miseq_test_batch_cluster/rep_set.fasta -uc miseq_test_batch_cluster/cluster_results.uc -id 0.97 -strand plus
```

The input to this function is the original list of dereplicated and ITS-extracted sequences, which is compared to the database of OTU representative sequences generated in the previous command. Essentially, in this step the original sequences are re-clustered by comparing to the representative sequences. Ideally we would run this function with the full list of sequences before dereplication, ITS extraction, and quality filtering (see http://drive5.com/usearch/manual/mapreadstootus.html for an explanation), but the need for ITS extraction precludes this approach for two reasons. First, and most importantly, running ITSx on reads that have not been quality filtered will result in many reads failing the homology search due to poor quality scores (leaving us discarding reads for unknown reasons as ITSx may also discard good reads with poor matches). Second, running ITSx on the full dataset takes considerably longer and may be untenable for very large datasets (without access to high computing power). The reads before ITS extraction may not map well (or at all) to the cluster representative sequences as they contain sequences of the 5.8S and LSU.

Step 8. Build OTU map
We can build an OTU map using the uc files generated during the sequence mapping step and the dereplication step using a perl script.

```
perl repo/fungal_ITS2_metabarcoding/combine_derep_and_clustering_uc_files.pl miseq_test_batch_cluster/cluster_results.uc miseq_test_batch_trim_merged_usearchFilter/seqs_ge_150.derep.uc miseq_test_batch_cluster/otus.txt
```

The script uses the two uc files to create an OTU map by searching for sequence ID’s that were clustered together. The output is written in the format of the mothur (Schloss et al. 2009) OTU map file, which was also adopted by QIIME (Caporasao et al. 2010) as a standard OTU map format. At this point we can remove the “barcodelabel=” portion of the sequence identifiers in the rep_set.fasta file and the OTU map file to be consistent with QIIME/mother format.

```
#Rewrite fasta and OTU file to QIIME style sequence ID format
perl repo/fungal_ITS2_metabarcoding/writeQiimeStyleRepSet.pl miseq_test_batch_cluster/rep_set.fasta miseq_test_batch_cluster/otus.txt miseq_test_batch_cluster/rep_set.qiimeStyle.fasta

perl repo/fungal_ITS2_metabarcoding/writeQiimeStyleOtuFile.pl miseq_test_batch_cluster/otus.txt miseq_test_batch_cluster/otus.qiimeStyle.txt
```

We’ll count the number of sequences and OTUs to make sure they match.

```
#Count
grep ">" miseq_test_batch_cluster/rep_set.qiimeStyle.fasta | wc
wc -l miseq_test_batch_cluster/otus.qiimeStyle.txt
```

Step 9. Assign taxonomy and make OTU table

The most common way to organize OTU by sample data for downstream analyses is in an “OTU table,” which is a table of sequence counts for each OTU (in rows) in each sample (in columns). The table often includes additional data provided in the column following the last sample, like a taxonomy string for each OTU. We will use QIIME to assign taxonomy by comparing our representative sequences to the UNITE reference database (Kõljalg et al. 2013) using blast. 
Please note that the UNITE database is a curated database of fungal ITS sequences. This means that most of the sequences in UNITE should be robust in terms of their identification and accuracy. BUT, the database may be missing species of fungi that are important in a specific dataset. AND, the database does not include anything other than fungi. If there is a possibility the primers used picked up something besides fungi (they probably did) then the sequences should also be compared to a comprehensive database like a recent version of the NCBI non-redundant nucleotide database, and positive matches to things like plants and animals should be removed. I usually accomplish this using blast, the program MEGAN, and perl scripts, and can provide that workflow if necessary.

```
#Assign taxonomy using QIIME
#or follow QIIME documentation for taxonomy assignment
qiime assign_taxonomy.py -i miseq_test_batch_cluster/rep_set.qiimeStyle.fasta -r sh_qiime_release_09/sh_refs_qiime_ver6_dynamic_09.02.2014.fasta -t sh_qiime_release_09/sh_taxonomy_qiime_ver6_dynamic_09.02.2014.txt -m blast -o miseq_test_batch_qiime
ls miseq_test_batch_qiime/
head miseq_test_batch_qiime/rep_set.qiimeStyle_tax_assignments.txt
```

Then build an OTU table.

``
#Make otu table
qiime make_otu_table.py -i miseq_test_batch_cluster/otus.qiimeStyle.txt -t miseq_test_batch_qiime/rep_set.qiimeStyle_tax_assignments.txt -o miseq_test_batch_qiime/otu_table.biom
```

Note that the file ending on the table is “.biom.” Biome (or biom) file format is a non-human readable file format that condenses metagenomic type data and sample metadata into a small file size. This is the standard table format for working in QIIME. To convert the biom table to a human readable format use convert_biom.py.

```
#Convert table from biom format to tab-delimited text format for viewing
#Check biom documentation for most recent command to accomplish this

biom convert -i miseq_test_batch_qiime/otu_table.biom -o miseq_test_batch_qiime/otu_table.txt -b --header-key taxonomy

#Look
less miseq_test_batch_qiime/otu_table.txt
```

The sequence data is now ready for downstream processing! Probably the most straight forward way to proceed is using QIIME, which has excellent documentation and many convenient scripts for subsampling sequences, calculating alpha-diversity and beta-diversity metrics, summing OTU sequence counts or relative abundance at various taxonomic levels, and some basic statistics functions. The ‘vegan’ package (Oksanen et al. 2015) in R is one of many R packages that also has some excellent analysis tools, but note that for most R functions the OTU table must be transposed.

References:
Caporaso JG, et al. (2010). QIIME allows analysis of high-throughput community sequencing data. Nature Methods 7(5): 335-336.

Edgar RC, (2013). UPARSE: highly accurate OTU sequences from microbial amplicon reads. Nature Methods. 10(10):996-998.

Edgar & Flyvbjerg (2015) Error filtering, pair assembly and error correction for next-generation sequencing. Bioinformatics. doi: 10.1093/bioinformatics/btv401

Kõljalg U, et al. (2013). Towards a unified paradigm for sequence-based identification of fungi. Molecular Ecology 22(21):5271-5277.

Oksanen J, et al. (2015). Vegan: community ecology package. http://CRAN.R-project.org/package=vegan

Schloss PD, et al. (2009). Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol 75(23):7537-41.

