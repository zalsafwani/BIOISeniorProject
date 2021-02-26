directory <- "/Users/zahraaalsafwani/Documents/htseq"

library("DESeq2")

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

#Pre-filtering the dataset
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

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

dds <- DESeq(dds)
(res <- results(dds))
summary(res)


mrna_res_df <- as.data.frame(res)
result_final_pos <- subset(mrna_res_df, res$log2FoldChange>0)
print(head(result_final_pos))

result_final_neg <- subset(mrna_res_df, res$log2FoldChange<0)
print(head(result_final_neg))

#write.csv(result_final_pos,file="DE_result_up.csv", quote=F)
#write.csv(result_final_neg,file="DE_result_down.csv", quote=F)

#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")

genesName <- c(row.names(res))
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=genesName,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

sessionInfo()
