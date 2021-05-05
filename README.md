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

### Re-run and make similar result to the GO anlysis
- After running the R script a file called DE_genes_network.txt, use this file to run GO in FunSet [FunSet website](http://funset.uno/). Then uplode the file and run, once the analysis is done download the csv result.
- Use it in R by replacing line 200.
- After running the code in R script for the GO analysis section, take the resulted file and copy the text into REVIGO [REVIGO website](http://revigo.irb.hr/).
- Change the orgnism to Homo sapiens, then run the analysis.

### Re-run and make similar result to the mRNA-miRNA network anlysis
- After running the R script a file called combine_ppi_miRNA_target_genes.csv
- Open Cytoscape File -> Import -> Network from File then choose the combine_ppi_miRNA_target_genes.csv
- Make column 1 as Not Imported, column 2 as Source Node, column 3 as Target Node, column 4 as Interaction Type, and column 5 as Source Node Attribute.
- Once the network created, File -> Import -> Table from File mRNA_regulatury.csv and make column 1 as Not Imported and column 3 (geneName) as a Key.
- File -> Import -> Table from File miRNA_regulatury.csv
- Choose the Style Node tab, change the Fill Color by selecting the Column to log2FoldChange and Mapping Type to Continuous Mapping.
- Change the Shape by selecting the Column to type and Mapping Type to Discrete Mapping and for the mRNA choose Ellipse and for miRNA choose Rectangle.
- Choose the Style Edge tab, change Stroke Color (Unselected) Column to interaction and Mapping Type to Discrete Mapping and for the PPI pick a green and for the mRNA-miRNA choose black. 
- Then choose Layout Group Attributies Layout -> Regulatury, then delate all of the nodes in the network that don't have Regulatury (blue nodes).
- To show the nodes info View -> Show Graphical Details.
- If we want to look at nodes with degree of connectivity > 100, Select tab then add Degree filter
- If we want to choose a node and look at the connected nodes, select the node from the Table Panel by right click then Select nodes from selectd rows. Select -> Nodes -> First Neighbors of Selected Nodes -> Undirected. Select -> Hid Unselected Nodes. Select -> Edges -> Show All Edges. To go back to the orginal network Select -> Show All Nodes and Edges

### Re-run and make similar result to the Somatic mutations anlysis
- After running the R script a file called genes_somatic_mutation.txt copy those genes. 
- Go to the cBioPortal [cBioPortal website](https://www.cbioportal.org/), filter the studies to only the Colorectal Adenocarcinoma then hit Query By Gene. Then past the genes list into the box and submit the query. 
- Click the Mutations tab to display the results.

