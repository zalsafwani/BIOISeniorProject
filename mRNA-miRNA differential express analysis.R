###########################################
# mRNA differential expressed analysis 
###########################################

# Set the directory, install and load the necessary library.
directory <- "~/Documents/htseq"
#BiocManager::install(c("DESeq2", "pathview"))
library(DESeq2)
library(pathview)

# Generate the HTSeqCount (mRNA) files as a DESeq dataset 
sampleFiles <- grep("counts",list.files(directory),value=TRUE)
sampleCondition <- sub("p[1,2,3,4](.*).htseq.counts","\\1",sampleFiles)
sampleName <- sub("(p[1,2,3,4].*).htseq.counts","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design= ~ condition)

# Pre-filtering the dataset
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)

# Run the differential expressed analysis for mRNA.
dds <- DESeq(dds)
(res <- results(dds))
summary(res)

# Save the result as data frame then filter it to the up_regulatory genes (mRNA)
mrna_res_df <- as.data.frame(res)
mRNA_res_filterd <- subset(mrna_res_df, mrna_res_df$padj<0.05)

# Map the gene symbol to a new column in the mRNA_res_filterd
setwd("~/Documents")
map <- read.table("Homo_sapiens.GRCh38.101.ensg2gene_name.txt",sep=" ",row.names=1)
geneName <- map[gsub("\\..*","",rownames(mRNA_res_filterd)),1]
mRNA_res_filterd <- cbind(mRNA_res_filterd, geneName)

# Add column based on DE mRNA UP or DOWN, type (mRNA), and Edge_type (mRNA-miRNA) 
# Last two column will be used when creating the network in Cytoscape
mRNA_res_filterd$Regulatury <- ifelse(mRNA_res_filterd$log2FoldChange>0, "UP", "DOWN")
mRNA_res_filterd$type <- "mRNA"
mRNA_regulatury <- subset(mRNA_res_filterd, select = c(0,2,6,7,8,9))
mRNA_regulatury$Edge_type <- "mRNA-miRNA"
write.csv(mRNA_regulatury,file="mRNA_regulatury.csv", quote=F)


###########################################
# KEEG for differential expressed mRNA
###########################################

# Formate them RNA_res_filterd in a specific way to run EGG pathway
mRNA_res_filterd_up <- subset(mRNA_res_filterd, mRNA_res_filterd$Regulatury == "UP")
upsymbol2eg <- id2eg(as.character(mRNA_res_filterd_up[,7]),category="symbol",org="Hs")
entrez_up <- upsymbol2eg[,2]
entrez_up <- as.numeric(unique(entrez_up[!is.na(entrez_up)]))

# Visualization of pathway and DEGs in KEGG pathway
# by assigning entrez ids to pvalues and remove entries with no entrez ids
p.values <- mRNA_res_filterd_up$pvalue
names(p.values) <-upsymbol2eg[,2]
p.values <- p.values[!is.na(names(p.values))]
pv.out <- pathview(gene.data = -log10(p.values), pathway.id = "05210", species = "hsa", out.suffix = "kegg_pathway")

# Formate them RNA_res_filterd in a specific way to run EGG pathway
mRNA_res_filterd_down <- subset(mRNA_res_filterd, mRNA_res_filterd$Regulatury == "DOWN")
upsymbol2eg <- id2eg(as.character(mRNA_res_filterd_down[,7]),category="symbol",org="Hs")
entrez_up <- upsymbol2eg[,2]
entrez_up <- as.numeric(unique(entrez_up[!is.na(entrez_up)]))

# Visualization of pathway and DEGs in KEGG pathway
# by assigning entrez ids to pvalues and remove entries with no entrez ids
p.values <- mRNA_res_filterd_down$pvalue
names(p.values) <-upsymbol2eg[,2]
p.values <- p.values[!is.na(names(p.values))]
pv.out <- pathview(gene.data = -log10(p.values), pathway.id = "05210", species = "hsa", out.suffix = "kegg_pathway_down")


###########################################
# miRNA differential expressed analysis 
###########################################

# Get the miRNA matrix and make it as DESeq dataset 
countdata <- read.table("miRNA_count_matrix.txt", header = T, row.names = 1)
countdata <- as.matrix(countdata)
(condition <- factor(c(rep("normal", 4), rep("tumor", 4))))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds_miRNA <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

#Pre-filtering the dataset
nrow(dds_miRNA)
dds_miRNA <- dds_miRNA[ rowSums(counts(dds_miRNA)) > 1, ]
nrow(dds_miRNA)
rld <- rlog(dds_miRNA, blind=FALSE)
head(assay(rld), 3)
par( mfrow = c( 1, 2 ) )
dds_miRNA <- estimateSizeFactors(dds_miRNA)

# Run the differential expressed analysis for miRNA.
dds_miRNA <- DESeq(dds_miRNA)
(res_miRNA <- results(dds_miRNA))
summary(dds_miRNA)

# Save the result as data frame then filter it to the down_regulatory genes (miRNA)
mirna_res_df <- as.data.frame(res_miRNA)
miRNA_res_filterd <- subset(mirna_res_df, mirna_res_df$padj<0.05)

# Add column based on DE miRNA UP or DOWN, type (miRNA), and Edge_type (mRNA-miRNA) 
# Last two column will be used when creating the network in Cytoscape
miRNA_res_filterd$Regulatury <- ifelse(miRNA_res_filterd$log2FoldChange>0, "UP", "DOWN")
miRNA_res_filterd$type <- "miRNA"
miRNA_regulatury <- subset(miRNA_res_filterd, select = c(0,2,6,7,8))
miRNA_regulatury$Edge_type <- "mRNA-miRNA"
write.csv(miRNA_regulatury,file="miRNA_regulatury.csv", quote = F)

# Get the table with miRNA to targeted gene for Homo sapiens and only keep miRNA column and the targeted gene
miRNAs_list <- read.csv("miRTarBase_SE_WR.csv")
miRNAs_list_Human <- subset(miRNAs_list, miRNAs_list$Species..miRNA. == "Homo sapiens")
miRNAs_targetGenes <- subset(miRNAs_list_Human, select = c(2,4))

# Remove duplicate rows
miRNAs_targetGenes <- miRNAs_targetGenes[!duplicated(miRNAs_targetGenes), ] 

# Fix the miRNA format and combine miRNA with targeted genes
targtedGenes_miRNA <- gsub("hsa-miR(.*)","hsa-mir\\1",miRNAs_targetGenes$miRNA)
mirna <- gsub("hsa(.*)-[3,5]p","hsa\\1",targtedGenes_miRNA)
miRNA <- cbind(miRNAs_targetGenes, mirna)
miRNA <- miRNA[!duplicated(miRNA), ] 
mirnas <- subset(miRNA, select = c(3,2))
mirnas$Edge_type <- "mRNA-miRNA"
write.csv(mirnas, file="miRNAs_targetGenes.csv", quote=F)


###########################################
# PPI file 
###########################################

# The file too large so I couldn't upload to githup
# Only select the column Official.Symbol.Interactor.A, Official.Symbol.Interactor.B
#ppi_data <- read.table("BIOGRID-ORGANISM-Homo_sapiens-4.3.196.tab3.txt", sep="\t", header = T, fill = TRUE)
#ppi <- subset(ppi_data, select = c(8,9))
#ppi$Edge_type <- "PPI"
#write.csv(ppi, ppi, quote = F)
ppi <- read.csv("ppi.csv")

# Add the genes to PPI if the gene symbol was in DE mRNA or DE miRNA based on their targeted genes
DE_miRNA <- data.frame(rownames(miRNA_regulatury))
PPI <- subset(ppi, ppi$Official.Symbol.Interactor.A %in% mRNA_regulatury$geneName | ppi$Official.Symbol.Interactor.B %in% mRNA_regulatury$geneName)
PPI <- subset(PPI, PPI$Official.Symbol.Interactor.A %in% mirnas$Target.Gene & mirnas$mirna %in% DE_miRNA$rownames.miRNA_regulatury.)

# Save the data as data frame then combine PPI with miRNA target genes into one file 
# It will be used to create the network in Cytoscape
miRNAs <- data.frame(unname(mirnas))
miRNAs$type <- "miRNA"
PPI <- data.frame(unname(PPI))
PPI$type <- "mRNA"
combine_ppi_miRNA_target_genes <- rbind(PPI, miRNAs)
write.csv(combine_ppi_miRNA_target_genes,file="combine_ppi_miRNA_target_genes.csv", quote = F)


###########################################
# Network analysis 
###########################################

# After creating the network and keeping only DE mRNAs and miRNAs save the result as csv
# Filter the results so only nodes with Degree >= 20 will be included 
# Then sort from high - low Degree
network_results <- read.csv("networkResultFixed.csv")
high_dgreeNode <- subset(network_results, network_results$Degree >= 20)
info_high_dgreeNode <- subset(high_dgreeNode, select = c(2,5,6,8,10))
info_high_dgreeNode <- info_high_dgreeNode[order(-info_high_dgreeNode$Degree),]

# Save both DE mRNAs and miRNAs form the network
# DE_genes_network will be used to make GO
# DE_miRNA_network will be used to catgorized miRNAs
DE_genes_network <- subset(network_results, network_results$type == "mRNA", select = c(5))
DE_miRNA_network <- subset(network_results, network_results$type == "miRNA", select = c(5,6))
write.table(DE_genes_network, file = "DE_genes_network.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(DE_miRNA_network$name, file = "DE_miRNA_net.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

###########################################
# GO analysis 
###########################################

# After running GO with DE_genes_network in FunSet and download the results as csv file
# Divide the GO into the two cluster to make the visualization in REVIGO
GO <- read.csv("final_GO.csv")
GO_cluster1 <- subset(GO, GO$Cluster == 0, select = c(2,4))
write.table(GO_cluster1, file = "GO_cluster1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
GO_cluster2 <- subset(GO, GO$Cluster == 1, select = c(2,4))
write.table(GO_cluster2, file = "GO_cluster2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

###########################################
# somatic mutation
###########################################

# Select mRNA if they have degree >= 100 to show their somatic mutation in cBioPortal
genes_somatic_mutation <- subset(info_high_dgreeNode, info_high_dgreeNode$type == "mRNA" & info_high_dgreeNode$Degree >= 100, select = c(1,2,3))
write.table(genes_somatic_mutation, file = "genes_somatic_mutation.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

###########################################
# catgorized miRNAs 
###########################################

# Read the gene_result.csv file that has miRNAs with their chromosomal locations
gene_result <- read.delim("~/Documents/gene_result.txt")
NCBI_miRNA_list <- subset(gene_result, select = c(6,7,9,10,11))

# Get the up-regulatory miRNAs with their chromosomal locations
DE_miRNA_network_up <- subset(DE_miRNA_network, DE_miRNA_network$Regulatury == "UP")
MIR_up <- gsub("hsa-mir-(.*)","MIR\\1",DE_miRNA_network_up$name)
MIR_MIRLET_up <- gsub("hsa-let-(.*)","MIRLET\\1",MIR_up)
MIR_up <- data.frame(MIR_MIRLET_up)
MIR_up <- apply(MIR_up,2,toupper)
chromosom_miRNA_up <- subset(NCBI_miRNA_list, select = c(1,4,5), NCBI_miRNA_list$Symbol %in% MIR_up)
write.csv(chromosom_miRNA_up, "chromosom_miRNA_up.csv", quote = F)

# Get the down-regulatory miRNAs with their chromosomal locations
DE_miRNA_network_down <- subset(DE_miRNA_network, DE_miRNA_network$Regulatury == "DOWN")
MIR_down <- gsub("hsa-mir-(.*)","MIR\\1",DE_miRNA_network_down$name)
MIR_MIRLET_down <- gsub("hsa-let-(.*)","MIRLET\\1",MIR_down)
MIR_down <- data.frame(MIR_MIRLET_down)
MIR_down <- apply(MIR_down,2,toupper)
chromosom_miRNA_down <- subset(NCBI_miRNA_list, select = c(1,4,5), NCBI_miRNA_list$Symbol %in% MIR_down)
write.csv(chromosom_miRNA_down, "chromosom_miRNA_down.csv", quote = F)

sessionInfo()


