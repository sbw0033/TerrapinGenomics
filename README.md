# TerrapinGenomics
Project with Matthew Wolak and Tonia Schwartz at Auburn University

This project starts with pre-processed raw reads. The files containing these reads are re-named and merged with all of the other in the script "TerrapinGenomicsPrepWork"

The major steps in moving from these raw reads to variant sites are as follows
1) Align the raw reads to a reference genome using the program BWA
2) Processing the raw alignments using Picard and Samtools
3) Calling variants using GATK's HaplotypeCaller

The GenomicsPrepWork.sh script uses the following code to add the sample ID to all fastqs and merge all the fastqs into two fastqs (1 forward 1 reverse)

```
#Copy all raw reads into the TestFastq directory
cp -a /hosted/biosc/SchwartzLab/WGS/Terrapin/WGS_BGI_03-2022/F22FTSUSAT0034_TURdgknR/soapnuke/clean/ /scratch/sbw0033/FastqCopies

#Make a loop that does the following- 1) Enters the sample directory 2) Adds the sample directory as a prefix to all files in that directory 3) exits the directory 4) does it again for every sample directory
cd /scratch/sbw0033/FastqCopies/clean/
for sampleID in *; do cd "${sampleID}"; for file in *.fq.gz; do mv "${file}" "${sampleID}_${file}"; done; cd ../; done #This does all the stuff mentioned above, written this way because I can't get indents right

#Make a loop that concatenates all the FW and reverse reads within each directory
for sampleID in *; do cd "${sampleID}"; cat *_1.fq.gz >> "${sampleID}_FW_reads.fq.gz"; cat *_2.fq.gz >> "${sampleID}_RV_reads.fq.gz"; cd ../; done #This does all the stuff mentioned above, written this way because I can't get indents right
```

### 1) RAW READ ALIGNMENT
When using BWA with GATK in mind as the eventual end goal, we need to be sure to include information about the sequencing platform, read group, and sample ID when running BWA.

BWA mem manual (http://bio-bwa.sourceforge.net/bwa.shtml)

### 2) RAW READ PROCESSING

Once we have all the alignments we need to do the following things
a) Convert .sam to .bam and sort it 
b) Ditch alignments where quality score is below 30
c) Remove all secondary alignments
d) [Add read group information so GATK doesn't freak out](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-)
e) [Mark duplicates- (Also so GATK doesn't freak out)](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)

### 3) CALL VARIANTS

We are going to use GATK's HaplotypeCaller to accomplish this
The BWAarray.sh script


