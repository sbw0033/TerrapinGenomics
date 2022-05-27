###Prefix for terrapin genomics Github is https://github.com/sbw0033/TerrapinGenomics.git

########################################################################################################
#Script to make sure that the reads have actually been properly trimmed for this project
########################################################################################################

#Transfer the painted turtle refernce genome and annotation to Easley
scp –r PaintedTurtleAnnotation.gtf.gz sbw0033@easley.auburn.edu:~ #Transfer the annotation
scp –r PaintedTurtleGenome.fa.gz sbw0033@easley.auburn.edu:~ #Transfer the reference genome
scp –r PaintedTurtleChr7.fa.gz sbw0033@easley.auburn.edu:/home/sbw0033/ #Transfer the reference genome

#Make a new scratch directory to test fastq stuff
cd /scratch/sbw0033
mkdir FastqCopiesTest
cd /scratch/sbw0033/FastqCopiesTest

#Enter the directory where all of the raw reads are located
cd /hosted/biosc/SchwartzLab/WGS/Terrapin/WGS_BGI_03-2022/F22FTSUSAT0034_TURdgknR/soapnuke/clean/

#Copy all raw reads into the TestFastq directory
cp -a /hosted/biosc/SchwartzLab/WGS/Terrapin/WGS_BGI_03-2022/F22FTSUSAT0034_TURdgknR/soapnuke/clean/ /scratch/sbw0033/FastqCopies

#Make a loop that does the following- 1) Enters the sample directory 2) Adds the sample directory as a prefix to all files in that directory 3) exits the directory 4) does it again for every sample directory
cd /scratch/sbw0033/FastqCopies/clean/
for sampleID in *; do cd "${sampleID}"; for file in *.fq.gz; do mv "${file}" "${sampleID}_${file}"; done; cd ../; done #This does all the stuff mentioned above, written this way because I can't get indents right

#Make a loop that concatenates all the FW and reverse reads within each directory
for sampleID in *; do cd "${sampleID}"; cat *_1.fq.gz >> "${sampleID}_FW_reads.fq.gz"; cat *_2.fq.gz >> "${sampleID}_RV_reads.fq.gz"; cd ../; done #This does all the stuff mentioned above, written this way because I can't get indents right

#These are all DNB-seq reads. Samples 2, 14, 76, and 93 were all sequenced at 40x depth instead of 20
