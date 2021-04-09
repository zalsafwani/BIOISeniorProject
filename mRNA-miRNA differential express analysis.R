# Set the directory, install and load the necessary library.
directory <- "/Users/zahraaalsafwani/Documents/htseq"
BiocManager::install(c("DESeq2", "pathview"))
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

# Shown are scatterplots using the log2 transform of normalized counts (left side)
#and using the rlog (right side)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)
sampleDists <- dist( t( assay(rld) ) )
sampleDists

# Run the differential expressed analysis for mRNA.
dds <- DESeq(dds)
(res <- results(dds))
summary(res)

# Save the result as data frame then filter it to the up_regulatory genes (mRNA)
mrna_res_df <- as.data.frame(res)
mRNA_res_filterd <- subset(mrna_res_df, mrna_res_df$padj<0.05)

# Map the gene symbol to a new column in the mRNA_res_filterd 
map <- read.table("Homo_sapiens.GRCh38.101.ensg2gene_name.txt",sep=" ",row.names=1)
geneName <- map[gsub("\\..*","",rownames(mRNA_res_filterd)),1]
mRNA_res_filterd <- cbind(mRNA_res_filterd, geneName)
write.csv(mRNA_res_filterd$geneName,file="mRNA_DE_result_up.csv", quote=F)

# Convert the 7th column to gene symbol
upsymbol2eg <- id2eg(as.character(mRNA_res_filterd[,7]),category="symbol",org="Hs")

# Save the second column of upsymbol2eg map into a vector called entrez_up
entrez_up <- upsymbol2eg[,2]
entrez_up <- as.numeric(unique(entrez_up[!is.na(entrez_up)]))

# Visualization of pathway and DEGs in KEGG pathway
# by assigning entrez ids to pvalues and remove entries with no entrez ids
p.values <- mRNA_res_filterd$pvalue
names(p.values) <-upsymbol2eg[,2]
p.values <- p.values[!is.na(names(p.values))]
pv.out <- pathview(gene.data = -log10(p.values), pathway.id = "05210", species = "hsa", out.suffix = "kegg_pathway")

###################################################################################################################################

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
miRNA_res_filterd <- subset(mirna_res_df, mirna_res_df$padj>0.05)
write.csv(miRNA_res_filterd,file="mirna_filter_res.csv", quote=F)

# Get the table with miRNA to targeted gene for Homo sapiens and only keep miRNA column and the targeted gene
miRNAs_list <- read.csv("miRTarBase_SE_WR.csv",)
miRNAs_list_Human <- subset(miRNAs_list, miRNAs_list$Species..miRNA. == "Homo sapiens")
miRNAs_targetGenes <- subset(miRNAs_list_Human, select = c(2,4))

#Remove duplicate rows
miRNAs_targetGenes <- miRNAs_targetGenes[!duplicated(miRNAs_targetGenes), ] 
row.names(miRNAs_targetGenes) < FALSE

targtedGenes_miRNA <- gsub("(.*)-[3,5]p","\\1",miRNAs_targetGenes$miRNA)
miRNA <- cbind(miRNAs_targetGenes, targtedGenes_miRNA)
miRNA <- miRNA[!duplicated(miRNA), ] 



sessionInfo()


