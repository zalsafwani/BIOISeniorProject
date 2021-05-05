# BIOISeniorProject
MicroRNAs and mRNAs differentially expressed analysis 

## Installation Language and Requirements:
- R as a language.
- Download the following packages and libraries: packages("BiocManager"), library(DESeq2), and library(pathview).

## Usage:
### executes the program
1. Downlowd the folowing files:
- htseq folder
- miRNA folder
- mRNA-miRNA differential express analysis.R
- Homo_sapiens.GRCh38.101.ensg2gene_name.txt
- miRNA_count_matrix.txt
- miRTarBase_SE_WR.csv
- ppi.csv
- networkResultFixed.csv
- final_GO.csv
- gene_result.txt
2. Open the mRNA-miRNA differential express analysis.R and set your working directory to the place where the files were downloded.
- Install "DESeq2" and "pathview" by uncommenting (line 9)
- Fix the directory to the htseq folder (line 6)
- Fix the directory to the folder that has all the files in step 1 (line 44)
4. Run the mRNA-miRNA differential express analysis.R
