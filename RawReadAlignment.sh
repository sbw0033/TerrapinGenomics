########################################################################################################
#Script to align raw reads for each sample to the painted turtle reference genome
########################################################################################################

#Load modules using specific version to ensure future compatibility
module load bwa/0.7.17

annotation="/home/sbw0033/TerrapinGenomics/Data/PaintedTurtleResources/PaintedTurtleAnnotation.gtf.gz"
reference="/home/sbw0033/TerrapinGenomics/Data/PaintedTurtleResources/PaintedTurtleGenome.fa.gz"
referenceChr7="/home/sbw0033/TerrapinGenomics/Data/PaintedTurtleResources/PaintedTurtleChr7.fa.gz"
RawReadsDir= "/scratch/sbw0033/FastqCopies/clean"

########################################################################################################
#First, need to build an index for the reference genome (painted turtle chr 7)
########################################################################################################
cd home/sbw0033/TerrapinGenomics/Data/PaintedTurtleResources
bwa index -a rb2 -p Chr7 PaintedTurtleChr7.fa
cd home/sbw0033/TerrapinGenomics/Data/PaintedTurtleResources/BWA_index_chr7 #Indexed files now live at /home/sbw0033/TerrapinGenomics/Data/PaintedTurtleResources/BWA_index_chr7

########################################################################################################
#Next, need to align the raw reads to the new reference index
########################################################################################################
#For BWA mem

for SAMPLE in 1..... ;#This should allow us to run bwa mem on all of our samples separately
do
	bwa mem \
	-t 48 \
	-R \

bwa mem -t 4 -R "@RG\tID:A8100\tSM:A8100" ../ref/chr18.fa A8100.chr18.R1.sickle.fastq A8100.chr18.R2.sickle.fastq > A8100.chr18.paired.sam

#Figure out where these annotations come from and how they were done
#



