#!/bin/bash

mkdir STAR-gencode

cd STAR-gencode

echo "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.primary_assembly.annotation.gtf.gz" >> url.txt
echo "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.primary_assembly.genome.fa.gz" >> url.txt

wget -i url.txt 

cd ..
