## 1. Installing STAR aligner with conda environment 


STAR version **2.7.6a** has been installed


### 1-1. Setting conda environment (environment.yml file)

check out conda docs for more info about environment management:    
https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html


```environment.yml
name: star
channels:
  - conda-forge
  - bioconda 
  - defaults 

dependencies: 
  - star =2.7.6a
  - r-base=4.0.2
  - r-tidyverse
  - r-data.table
  - r-ggplot2
  - r-markdown
  - r-pheatmap
  - bioconductor-deseq2
  - bioconductor-annotationhub
  - bioconductor-tximport
  - bioconductor-rsubread
  - bioconductor-apeglm

```

### 1-2. Installing the conda enviornment using environment.yml 


```terminal


conda env create -f environment.yml

```

### 1-3. Activating the conda environment 

```terminal
conda activate star

```

## 2. Downloading reference files 

- Reference genome sequences (FASTA files): https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ ---> hg19.fa.gz ---> unzip 
- Annotations (GTF file): https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/  ---> hg19.knownGene.gtf.gz ---> unzip 

GTF files have multiple versions such as **ncbiRefSeq, refGene, ensGene, and knownGene**. For more info, see below:  

For the human assembly hg19/GRCh37 and mouse mm9/NCBI37: What is the difference between UCSC Genes, the "GENCODE Gene Annotation" track and the "Ensembl Genes" track?    

The "UCSC Genes" track, also called "Known Genes", is available only on assemblies before hg38. It was built with a gene predictor developed at UCSC. This gene predictor uses protein, EST and cDNA annotations to derive a relatively restricted gene transcript set. The software is no longer in use and there are no plans to release the track on newer human assemblies. It was last used for the mm10 mouse assembly.    

The "GENCODE Gene Annotation" track contains data from all versions of GENCODE. "Ensembl Genes" track contains just a single Ensembl version. See the previous question for the differences between Ensembl and GENCODE.   

(https://genome.ucsc.edu/FAQ/FAQgenes.html#hg19)


**Note**: It has to be checked out whether the version of my STAR aligner is compatible with the reference genome



## 3. Running STAR (2-pass mode)

### Reference documents 
- PMID 23104886
- docs and manual: https://github.com/alexdobin/STAR


### 3-1. Generating genome indexing (STARindex.sh) 

- It has to be run in the directory where STARindex.sh, the genome FASTA file (e.g. hg19.fa), and the GTF file (e.g. hg19.knownGene.gtf) locate 
- Set **"--runMode genomeGenerate"** 
- The output index files are saved in ./genomegen

```STARindex.sh


#!bin/bash 


fasta="hg19.fa"
gtf="hg19.knownGene.gtf"

mkdir genomegen

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./genomegen --genomeFastaFiles $fasta --sjdbGTFfile $gtf 


```

Conducting genome indexing by executing the STARindex.sh file 


```terminal


(star) mira@mira-MS-7C90:~/Documents/programming/Bioinformatics/STAR-test/Hs_Hg19$ bash STA
Rindex.sh
Nov 01 22:42:18 ..... started STAR run
Nov 01 22:42:18 ... starting to generate Genome files
Nov 01 22:42:57 ..... processing annotations GTF
Nov 01 22:43:09 ... starting to sort Suffix Array. This may take a long time...
Nov 01 22:43:22 ... sorting Suffix Array chunks and saving them to disk...
Nov 01 22:59:33 ... loading chunks from disk, packing SA...
Nov 01 23:00:23 ... finished generating suffix array
Nov 01 23:00:23 ... generating Suffix Array index
Nov 01 23:03:30 ... completed Suffix Array index
Nov 01 23:03:30 ..... inserting junctions into the genome indices
Nov 01 23:05:35 ... writing Genome to disk ...
Nov 01 23:05:36 ... writing Suffix Array to disk ...
Nov 01 23:05:49 ... writing SAindex to disk
Nov 01 23:05:51 ..... finished successfully


```



### 3-2. Running alignment (star_test.sh)


#### Basic parameters
- **--runThreadN**: number of threads
- **--runMode**: genomeGenerate (for indexing) or alignReads (for alignment. by default)
- **--genomeDir**: path to genomeDir (contains index files)
- **--genomeFastaFiles**: path to ref genome fasta files (e.g. Homo_sapiens.GRCh38.dna.primary_assembly.fa)
- **--sjdbGTFfile**: path to annotations.gtf (reference transcripts)
- **--sjdbOverhang**: ReadLength-1 (or 100 is also descent)
- **--readFilesIn**: path/name of input files
- **--outFileNamePrefix**: path/name of output files

#### Advanced parameters
- **--outFilterType BySJout**: reduces number of spurious junctions
- **--outFilterMultimapNmax 20**: if more than this many multimappers, consider unmapped
- **--alignSJoverhangMin 8**: min overhang for unannotated junctions
- **--alignSJDBoverhangMin 1**: min overhang for annotated junctions
- **--outFilterMismatchNmax 999**: max mismatches per pair
- **--outFilterMismatchNoverReadLmax 0.04**: max mismatches per pair relative to read length
- **--alignIntronMin 20**: min intron length
- **--alignIntronMax 1000000**: max intron length
- **--outSAMunmapped None**: do not report aligned reads in output
- **--outSAMtype BAM SortedByCoordinate**: output sorted by coordinate 
- **--quantMode GeneCounts**: ouput count number per gene (ReadsPerGene.out.tab file)
- **--twopassMode Basic**: perform 2-pass mapping automatically  
- **--chimOutType Junctions**: detect chimeric alignments by giving output Chimeric.out.junction files


Note: 
- Once indexing is correctly completed, keeping the **--genomeFastaFiles** parameter causes an error when running alignment. Avoid using the parameter for alignment.
- It is recommended to set **ref_path** to an absolute path due to an error


Create star_test.sh file 

```star_test.bash
#!bin/bash


# Define paths file names depending on your current data

mkdir output 

input_path="../rawdata"
ref_path="/home/mira/Documents/programming/Bioinformatics/STAR-test/Hs_Hg19" # absolute path
genome_fasta="hg19.fa"
GTF="hg19.knownGene.gtf"

cd rawdata

input_files=$(ls)

cd ../output



for read in $input_files 
do 
    STAR --runThreadN 8 --runMode alignReads --genomeDir ${ref_path}/genomegen --sjdbGTFfile ${ref_path}/${GTF} -sjdbOverhang 100 --readFilesIn ${input_path}/${read} --outFileNamePrefix ${read} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic --chimOutType Junctions

done


```

Execute the star_test.sh file 

```terminal

bash star_test.sh

```

Running


```

(star) mira@mira-MS-7C90:~/Documents/programming/Bioinformatics/STAR-test$ bash star_test.sh
Nov 01 23:35:27 ..... started STAR run
Nov 01 23:35:27 ..... loading genome
Nov 01 23:35:42 ..... processing annotations GTF
Nov 01 23:35:46 ..... inserting junctions into the genome indices
Nov 01 23:36:36 ..... started 1st pass mapping
Nov 01 23:39:01 ..... finished 1st pass mapping
Nov 01 23:39:01 ..... inserting junctions into the genome indices
Nov 01 23:40:06 ..... started mapping
Nov 01 23:42:53 ..... finished mapping
Nov 01 23:42:54 ..... started sorting BAM
Nov 01 23:43:17 ..... finished successfully
Nov 01 23:43:18 ..... started STAR run
Nov 01 23:43:18 ..... loading genome
Nov 01 23:43:34 ..... processing annotations GTF
Nov 01 23:43:39 ..... inserting junctions into the genome indices
Nov 01 23:44:30 ..... started 1st pass mapping
Nov 01 23:46:40 ..... finished 1st pass mapping
Nov 01 23:46:40 ..... inserting junctions into the genome indices
Nov 01 23:47:46 ..... started mapping
Nov 01 23:50:15 ..... finished mapping
Nov 01 23:50:16 ..... started sorting BAM
Nov 01 23:50:36 ..... finished successfully
Nov 01 23:50:37 ..... started STAR run
Nov 01 23:50:37 ..... loading genome
Nov 01 23:50:54 ..... processing annotations GTF
Nov 01 23:50:59 ..... inserting junctions into the genome indices
Nov 01 23:51:50 ..... started 1st pass mapping
Nov 01 23:53:39 ..... finished 1st pass mapping
Nov 01 23:53:39 ..... inserting junctions into the genome indices
Nov 01 23:54:43 ..... started mapping
Nov 01 23:56:49 ..... finished mapping
Nov 01 23:56:50 ..... started sorting BAM
Nov 01 23:57:07 ..... finished successfully
Nov 01 23:57:08 ..... started STAR run
Nov 01 23:57:08 ..... loading genome
Nov 01 23:57:24 ..... processing annotations GTF
Nov 01 23:57:28 ..... inserting junctions into the genome indices
Nov 01 23:58:20 ..... started 1st pass mapping
Nov 02 00:01:38 ..... finished 1st pass mapping
Nov 02 00:01:38 ..... inserting junctions into the genome indices
Nov 02 00:02:44 ..... started mapping
Nov 02 00:06:28 ..... finished mapping
Nov 02 00:06:29 ..... started sorting BAM
Nov 02 00:06:55 ..... finished successfully
Nov 02 00:06:57 ..... started STAR run
Nov 02 00:06:57 ..... loading genome
Nov 02 00:07:14 ..... processing annotations GTF
Nov 02 00:07:18 ..... inserting junctions into the genome indices
Nov 02 00:08:09 ..... started 1st pass mapping
Nov 02 00:09:46 ..... finished 1st pass mapping
Nov 02 00:09:47 ..... inserting junctions into the genome indices
Nov 02 00:10:48 ..... started mapping
Nov 02 00:12:39 ..... finished mapping
Nov 02 00:12:40 ..... started sorting BAM
Nov 02 00:12:55 ..... finished successfully
Nov 02 00:12:56 ..... started STAR run
Nov 02 00:12:56 ..... loading genome
Nov 02 00:13:13 ..... processing annotations GTF
Nov 02 00:13:17 ..... inserting junctions into the genome indices
Nov 02 00:14:06 ..... started 1st pass mapping
Nov 02 00:16:30 ..... finished 1st pass mapping
Nov 02 00:16:30 ..... inserting junctions into the genome indices
Nov 02 00:17:31 ..... started mapping
Nov 02 00:20:14 ..... finished mapping
Nov 02 00:20:15 ..... started sorting BAM
Nov 02 00:20:36 ..... finished successfully

```

## 4. Extracting transcript read counts from BAM files with featureCounts in R (Rsubread package)

- Annotation was done with Ensembl Hg19 GTF (downloadable on: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes)
- The featureCounts is available via Unix and R packages. In this workflow, R package was used.  
- Check out featureCounts docs (http://subread.sourceforge.net) for more info




```r

# Loading packages 
library(tidyverse)
library(Rsubread)
library(AnnotationHub)

# AnnotationHub setup 
AnnotationSpecies <- "Homo sapiens"  # Assign your species 
MyCache <- "../rawdata/AnnoCache"
ah <- AnnotationHub(hub=getAnnotationHubOption("URL"))   # Bring annotation DB


# AnnotationHub run

# Explore your OrgDb object with following accessors:
# columns(AnnpDb)
# keytypes(AnnoDb)
# keys(AnnoDb, keytype=..)
# select(AnnoDb, keys=.., columns=.., keytype=...)
AnnoKey <- keys(AnnoDb, keytype="ENSEMBLTRANS")
# Note: Annotation has to be done with not genome but transcripts 
AnnoDb <- select(AnnoDb, 
                 AnnoKey,
                 keytype="ENSEMBLTRANS",
                 columns="SYMBOL")
head(AnnoDb)




# Defining featureCounts parameters (has to be done by users)
Samples <- c("Mock_72hpi_S1",
             "Mock_72hpi_S2",
             "Mock_72hpi_S3",
             "SARS-CoV-2_72hpi_S7",
             "SARS-CoV-2_72hpi_S8",
             "SARS-CoV-2_72hpi_S9")


NameTail=".fastq.gzAligned.sortedByCoord.out.bam"


BAMInputs <- c()   # Path to Input files 

for (i in 1:length(Samples)) {
    BAMInputs[i] <- paste0("./output/", 
                           Samples[i],
                           NameTail)
}



# "mm10", "mm9", "hg38", or "hg19"
annot.inbuilt="hg19" 

# annotation data file such as a data frame or a GTF file
annot.ext="./featurecountdata/hg19.ensGene.gtf"

# What to annotate 
GTF.attrType="transcript_id"




# Define a function extracting transcript counts from a BAM file ans saving as a data frame 
ExtractCounts_fn <- function(BAMInput) {

    # Run featureCounts() function
    FC <- featureCounts(BAMInput,
                    annot.inbuilt=annot.inbuilt,
                    annot.ext=annot.ext,
                    GTF.attrType=GTF.attrType,
                    isGTFAnnotationFile=TRUE,
                    nthreads=nthreads,
                    verbose=TRUE)

    # Extract a count matrix
    FCmatrix <- FC$counts

    # Save as a data frame
    Transcripts=rownames(FCmatrix)
    Counts=FCmatrix[,1]
    CountTable <- data.frame(Transcripts, Counts)
    
    return(CountTable)
}


# Initialize CountTable data frame with the first BAM file
CountTable <- ExtractCounts_fn(BAMInputs[1])



# Combine count data from the second to the last BAM file to the CountTable data frame
for (i in 2:length(Samples)) {
         ct <- ExtractCounts_fn(BAMInputs[i])
         CountTable <- full_join(CountTable, ct, by="Transcripts")
}


# Remove undetected transcripts from the CounTable data frame
CountTable <- CountTable[rowSums(CountTable[, -1]) != 0,]

# Change column names 
ColRename <- c("Transcripts", Samples)
colnames(CountTable) <- ColRename

# Add gene name annotation
CountTable <- right_join(AnnoDb, 
                         CountTable, 
                         by=c("ENSEMBLTRANS"="Transcripts")) 
colnames(CountTable) <- c("Transcript", "Gene", Samples)





# Create a directory and save as a csv file
dir.create("./csv")
write.csv(CountTable, "./csv/read_count.csv") 
```

