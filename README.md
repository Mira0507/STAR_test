## STAR aligner 

### 1. Raw data
- Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157852
- Same raw data as the previous [Salmon workflow](https://github.com/Mira0507/salmon_test)

### 2. workflow_ens.md
- Workflow from raw data to read count 
- Alignment with **STAR** producing BAM files and counting reads with **featureCounts** 

### 3. featureCounts_ens.Rmd 
- Counting reads 
- Simple statistics 
- Downstream DE analysis
