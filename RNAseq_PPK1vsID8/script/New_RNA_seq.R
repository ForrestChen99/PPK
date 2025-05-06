setwd('/work/cwt/ID8/4_RNA_seq')
library(DESeq2)
library(EnhancedVolcano)

rna_count = read.delim("featureCounts_merged_count.annot.tsv")
rna_count[1:5,1:5]
rna_count1 = dplyr::select(rna_count,gene_name,ID8.1,ID8.2,ID8.3,ID8.PK.1,ID8.PK.2,ID8.PK.3)
rna_count2 = distinct(rna_count1,gene_name,.keep_all = T)
rownames(rna_count2) = rna_count2$gene_name
rna_count2 = rna_count2[-1]

phe = data.frame(group = c(rep('ID8',3),rep('ID8_PK',3)),
                 sample = c('ID8.1','ID8.2','ID8.3','ID8.PK.1','ID8.PK.2','ID8.PK.3'))
rownames(phe) = phe$sample
phe$group = factor(phe$group,levels = c('ID8','ID8_PK'))


dds = DESeqDataSetFromMatrix(countData = rna_count2,
                             colData = phe,
                             design = ~group,
                             tidy = F)
head(dds)
dds = dds[rowSums(counts(dds))>=10,]
dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=F))


EnhancedVolcano(res,
                lab = rownames(res),
                labSize = 0,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 2,
                ylab = bquote(~-Log[2]~ 'fold change'))

res1<-as.data.frame(res)

resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

plotCounts(dds, gene=which.min(res$padj), intgroup="group")

res2<- results(dds,contrast=c('group','ID8_PK','ID8'), altHypothesis="greaterAbs",lfcThreshold = 1)

res2[which(res2$padj<0.05),]
res3 = as.data.frame(res2)

