setwd('/work/cwt/ID8/4_RNA_seq')
library(DESeq2)
library(EnhancedVolcano)

##DESEQ2
{
  deseqdt = readRDS('DESeq2.rds')
  
  list.files()
  rna_count = read.delim("featureCounts_merged_count.annot.tsv")
  
  rna_count1 = dplyr::select(rna_count,gene_name,ID8.1,ID8.2,ID8.3,ID8.PK.1,ID8.PK.2,ID8.PK.3)
  #rownames(rna_count1) = rna_count1$gene_name
  
  rna_count2 = distinct(rna_count1,gene_name,.keep_all = T)
  rownames(rna_count2) = rna_count2$gene_name
  rna_count2 = rna_count2[-1]
  
  #批次效应检验
  dat = as.data.frame(log2(edgeR::cpm(rna_count2)+1))
  color <- brewer.pal(7, 'RdGy')
  boxplot(dat, col=color, las=2)
  
  
  phe = data.frame(group = c(rep('ID8',3),rep('ID8_PK',3)),
                   sample = c('ID8.1','ID8.2','ID8.3','ID8.PK.1','ID8.PK.2','ID8.PK.3'))
  rownames(phe) = phe$sample
  phe$group = factor(phe$group,levels = c('ID8','ID8_PK'))
  
  
  #DEseq2
  dds = DESeqDataSetFromMatrix(countData = rna_count2,
                               colData = phe,
                               design = ~group,
                               tidy = F)
  dds = dds[rowSums(counts(dds))>=10,]
  dds
  dds <- DESeq(dds)
  
  res <- results(dds)
  head(results(dds, tidy=F))
  summary(res)
  
  
  resordered <- res[order(res$padj),]
  head(resordered)
  sum(resordered$padj < 0.1, na.rm=TRUE)
  
  
  res = na.omit(res)
  res[res$padj == 0] = 0.001
  
  EnhancedVolcano(res,
                  lab = rownames(res),
                  labSize = 0,
                  x = 'log2FoldChange',
                  y = 'padj',
                  pCutoff = 0.05,
                  FCcutoff = 2,
                  ylab = bquote(~-Log[2]~ 'fold change'))
  
  colnames(res)
  
  {
    res1<-as.data.frame(res)
    res1= na.omit(res1)
    res1 = dplyr::filter(res1,padj != 0)
    
    
    threshold<-as.factor(
      (res1$log2FoldChange>2.5|res1$log2FoldChange<(-2.5)) & 
        res1$padj<0.05)
    
    
    ggplot(res1,aes(x= log2FoldChange,
                    y= -1*log10(padj),colour=threshold),palet )+xlab("log2 fold-change")+ylab("-log10 Padj")+geom_point()+
      theme_classic()+scale_color_manual(values = c('#4DBBD57F','#E64B357F'))
    
    write.csv(res1,file = 'Deseq2result.csv')
  }
  
}

##ANALYSIS OF DOWN/UP REGULATED GENES
up = dplyr::filter(res1,log2FoldChange>0.585 & adj.P.Val < 0.05)
down = dplyr::filter(res1, log2FoldChange<(-0.585) & padj < 0.05)

ego_up = enrichGO(
  gene = rownames(up),
  keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
  OrgDb = org.Mm.eg.db,
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = FALSE
)
ego_down = enrichGO(
  gene = rownames(down),
  keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
  OrgDb = org.Mm.eg.db,
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = FALSE
)
dotplot(ego_down)
dotplot(ego_up)
plotdata_up = ego_up@result
plotdata_up$group = 'up'
plotdata_up = plotdata_up %>% dplyr::do(head(.,n=5)) 

plotdata_down = ego_down@result
plotdata_down$group = 'down'
plotdata_down = plotdata_down %>% dplyr::do(head(.,n=5)) 

plotdata = rbind(plotdata_up,plotdata_down)

plotdata = plotdata %>% mutate(padj1 = ifelse(group == 'up',-log10(p.adjust),log10(p.adjust))) %>% arrange(group,padj1)
plotdata$group = factor(plotdata$group,levels = c('up','down'))

ggplot(plotdata,aes(x = reorder(Description,padj1),y = padj1))+
  geom_bar(aes(fill=group),stat="identity")+
  scale_fill_manual(values = c('#91D1C2FF','#F39B7FFF'))+
  labs(x = '',y = '-log10(adjusted Pvalue)',fill = '')+
  theme_classic()+ggtitle('')+coord_flip()
ggsave(filename = 'RNAseq_goenrich.pdf',width=8,height = 3)


##投射gene feature回tumor featureplot
ppk1.gene = up %>% top_n(40,wt = log2FoldChange) %>% rownames()

rows_to_remove <- grepl("^Gm", rownames(down))
down <- down[!rows_to_remove, ]

id8wt.gene = down %>% top_n(-40,wt = log2FoldChange) %>% rownames()

DefaultAssay(tumor_ident) = 'RNA'
tumor_ident = AddModuleScore(object = tumor_ident,
                             features = as.list(ppk1.gene[ppk1.gene %in% rownames(tumor_ident)]),
                             name = 'ppk1.gene.0212')
tumor_ident = AddModuleScore(object = tumor_ident,
                             features = as.list(id8wt.gene[id8wt.gene %in% rownames(tumor_ident)]),
                             name = 'id8wt.gene1.0212')

FeaturePlot(tumor_ident,features = c('ppk1.gene.02121','id8wt.gene1.02121'),ncol=2,
            cols = c('grey','orange','red'),order = T,raster = T)
ggsave(filename = 'feature_rnaseq.pdf',width = 8,height=4)

tumor_ident
table(tumor_ident$group)
tnk
