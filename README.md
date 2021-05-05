# BIOISeniorProject
MicroRNAs and mRNAs differentially expressed analysis using DESeq2

The goal of this project is to run differential expression analysis for mRNAs and miRNAs data in Colorectal cancer (CRC) and see how miRNAs and their targeted genes could correlate with CRC. Patients with CRC have high mortality rates because they are diagnosed at an advanced stage. Identifying miRNAs–mRNAs that have a role in cellular pathways, such as apoptosis, could help identify markers for cancer, especially of CRC. Understanding expressions of miRNAs with their targeted mRNAs could help scientists find biomarkers that could help early diagnosis and better treatments of CRC with higher survival rates.

## Installation Language and Requirements:
- R
- Packages and libraries: packages("BiocManager"), library(DESeq2), and library(pathview)

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

### session info
More information about the session info are provided in session_info file

