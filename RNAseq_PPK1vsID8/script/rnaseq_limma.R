library(limma)

##expression data
rna_count = read.delim("featureCounts_merged_count.annot.tsv")
rna_count[1:5,1:5]
rna_count1 = dplyr::select(rna_count,gene_name,ID8.1,ID8.2,ID8.3,ID8.PK.1,ID8.PK.2,ID8.PK.3)
rna_count2 = distinct(rna_count1,gene_name,.keep_all = T)
rownames(rna_count2) = rna_count2$gene_name
rna_count2 = rna_count2[-1]

#removing of unwanted genes
rna_count3 = rna_count2[!grepl("^Rp[sl]",rownames(rna_count2),ignore.case = F),]
rna_count3 = rna_count3[!grepl("^mt-",rownames(rna_count3),ignore.case = F),]
rna_count3 = rna_count3[!grepl("^Hsp",rownames(rna_count3),ignore.case = F),]


#design matrix
group_list = c(rep('ID8',3),rep('ID8_PK',3)) %>% factor() %>% relevel(ref = "ID8")
design = model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(expr_normalize)


#edgeR: data processing
library(edgeR)
dge = DGEList(counts = rna_count3) #rna_count3 remove low quality gene
keep = filterByExpr(dge)
dge = dge[keep,,keep.lib.sizes=FALSE]

x = dge
group = as.factor(c(rep('ID8',3),rep('ID8_PK',3)))
x$samples$group = group

cpm = cpm(x)
lcpm = cpm(x,log=TRUE,prior.count = 2)

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#differentially expressed genes
design

fit = lmFit(lcpm,design)
fit = eBayes(fit,trend = T)#
topTable(fit,coef = 2)
dt1 = topTable(fit,adjust='BH',number = Inf)

dt1['Mal2',]

allDiff = topTable(fit,coef=2,adjust='fdr',number = Inf)
allDiff['Mal2',]
lcpm['Mal2',]
data = allDiff
write.csv(data,'DEG.csv')
data$significant="stable"
data$significant[data$logFC>=0.585 & data$adj.P.Val <0.05]="up"
data$significant[data$logFC<= -0.585 & data$adj.P.Val <0.05]="down"

ggplot(data,aes(logFC,-1*log10(adj.P.Val)))+
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#4CBBD5","grey","#f39b7f"))+
  geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-0.5,0.5),linetype=4,size=0.3)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2(Fold Change)",y="-log10(Adjusted Pvalue)")
ggsave(filename = 'Volcano.pdf',width = 6, height = 4)

#ggrepel gene name selection
if(F){
  select.lcpm <- (diff$AveExpr > 5 )
  table(select.lcpm)
  
  select.log2FC <- abs(diff$logFC) >0.585
  table(select.log2FC)
  select.qval <- (diff$adj.P.Val< 0.05)
  table(select.qval)
  
  select.vec=(select.lcpm & select.log2FC & select.qval) 
  table(select.vec)
  
  degs.list=as.character(diff$gene)[select.vec]
  
  label.deg=sample(degs.list,20)
  p = ggplot(data,aes(logFC,-1*log10(adj.P.Val)))+
    geom_point(aes(color=significant),size=0.8)+theme_classic()+
    scale_color_manual(values = c("#4CBBD5","grey","#f39b7f"))+
    geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
    geom_vline(xintercept = c(-0.5,0.5),linetype=4,size=0.3)+
    theme(title=element_text(size = 18),text = element_text(size=18))+
    labs(x="log2(Fold Change)",y="-log10(Adjusted Pvalue)")
  
  #data_selected <- data[label.deg,]
  #p+geom_label_repel(data=data_selected,aes(label=rownames(data_selected)))
  }

##画个热图

up.gene = data %>% top_n(10,wt = logFC) %>% rownames()
down.gene = data %>% top_n(-12,wt = logFC) %>% rownames()
down.gene = down.gene[-c(9,12)]

gene.select = c(down.gene,up.gene)

exprSet.map=lcpm[gene.select,]

annotation_col1 = data.frame(Group =as.factor(c(rep('ID8-wt',3),rep('ID8-PPK1',3))))
rownames(annotation_col1)=colnames(lcpm)
ann_colors = list(Group = c(`ID8-wt`="#91D1C27F", `ID8-PPK1`="#3C54887F"))
p = pheatmap::pheatmap(exprSet.map, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   annotation_col =annotation_col1,
                   annotation_colors = ann_colors,
                   show_colnames=T,
                   #scale = "row", #以行来标准化，这个功能很不错
                   color = rev(colorRampPalette(brewer.pal(2,'RdBu'))(50)))
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p,filename = 'heatmap.pdf')

#############Gene ontology analysis###############
{
  up = dplyr::filter(data,logFC>=0.585 & adj.P.Val < 0.05)
  down = dplyr::filter(data, logFC<=(-0.585) & adj.P.Val < 0.05)
  
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
  ego_up_all = enrichGO(
    gene = rownames(up),
    keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
    OrgDb = org.Mm.eg.db,
    ont = 'ALL',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = FALSE
  )
  ego_down_all = enrichGO(
    gene = rownames(down),
    keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
    OrgDb = org.Mm.eg.db,
    ont = 'ALL',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = FALSE
  )
  dotplot(ego_down_all,split = 'ONTOLOGY')
  dotplot(ego_up_all,split = 'ONTOLOGY')
  
  plotdata_up = ego_up_all@result
  plotdata_up$group = 'up'
  plotdata_up1 = plotdata_up %>% group_by(ONTOLOGY) %>% top_n(-5,wt = p.adjust)
  
  p1 = ggplot(plotdata_up1,aes(x = reorder(Description,p.adjust),y = Count))+
    geom_bar(aes(fill=p.adjust),stat="identity",position = "dodge")+
    facet_grid(ONTOLOGY~.,scales = "free",space = "free")+
    scale_fill_gradient(low = '#EF8A62', high = '#67A9CF')+
    labs(x = '',y = 'Count',fill = '')+
    theme_classic()+ggtitle('ID8-PPK1')+coord_flip()
  
  plotdata_down = ego_down_all@result
  plotdata_down$group = 'down'
  plotdata_down1 = plotdata_down %>% group_by(ONTOLOGY) %>% top_n(-5,wt = p.adjust)
  
  p2 = ggplot(plotdata_down1,aes(x = reorder(Description,p.adjust),y = Count))+
    geom_bar(aes(fill=p.adjust),stat="identity",position = "dodge")+
    facet_grid(ONTOLOGY~.,scales = "free",space = "free")+
    scale_fill_gradient(low = '#EF8A62', high = '#67A9CF')+
    labs(x = '',y = 'Count',fill = '')+
    theme_classic()+ggtitle('Downregulated')+coord_flip()
  p = p1 + p2
  ggsave(p,filename = 'GO.pdf',height = 5,width = 12)
  
  #plotdata = rbind(plotdata_up1,plotdata_down1)
  #plotdata = plotdata %>% mutate(padj1 = ifelse(group == 'up',-log10(p.adjust),log10(p.adjust))) %>% arrange(group,padj1)
  #plotdata$group = factor(plotdata$group,levels = c('up','down'))
}

###GSEA enrichment analysis
#install.packages('msigdbr')
library(msigdbr)
library(fgsea)


m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

colnames(data)
data$gene = rownames(data)
genelist = data %>% dplyr::filter(logFC>0) %>% arrange(desc(logFC)) %>%
  dplyr::select(gene,logFC)

genelist1 = deframe(genelist)

fgseaRes = fgsea(fgsea_sets,stats = genelist1)

ggplot(fgseaRes %>% filter(padj < 0.1) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES)) +
  coord_flip() +
  scale_fill_gradient(low = '#FDE0DD', high = '#C51B8A')+
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways NES from GSEA") +
  theme_minimal() ####以1.5进行绘图填色
ggsave(filename = 'KEGG_gsea.pdf',width= 10,height=5,path = '/work/cwt/ID8/4_RNA_seq/')
brewer.pal(2,'RdPu')

p1 = plotEnrichment(fgsea_sets[["KEGG_PATHWAYS_IN_CANCER"]],
               genelist1) + labs(title="KEGG_PATHWAYS_IN_CANCER")
p2 = plotEnrichment(fgsea_sets[["KEGG_CELL_ADHESION_MOLECULES_CAMS"]],
               genelist1) + labs(title="KEGG_CELL_ADHESION_MOLECULES_CAMS")
p = p1+p2
pdf('fgsea_kegg.pdf',width = 8,height = 3)
p
dev.off()

fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))


#############Reactome pathway###############
m_df = msigdbr(species = "Mus musculus")
m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat) %>% print(n=20)

m_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes = fgsea(fgsea_sets,stats = genelist1)

ggplot(fgseaRes %>% filter(padj < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES)) +
  coord_flip() +
  scale_fill_gradient(low = '#EFEDF5', high = '#756BB1')+
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Reactome pathways NES from GSEA") +
  theme_minimal() ####以1.5进行绘图填色

brewer.pal(2,'Purples')
ggsave(filename = 'Reactome pathways NES from GSEA.pdf',width = 10,height = 5)

p1 = plotEnrichment(fgsea_sets[["REACTOME_GPCR_LIGAND_BINDING"]],
                    genelist1) + labs(title="REACTOME_GPCR_LIGAND_BINDING")
p2 = plotEnrichment(fgsea_sets[["REACTOME_INTERFERON_GAMMA_SIGNALING"]],
                    genelist1) + labs(title="REACTOME_INTERFERON_GAMMA_SIGNALING")
p3 = plotEnrichment(fgsea_sets[["REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"]],
                    genelist1) + labs(title="REACTOME_INTERFERON_ALPHA_BETA_SIGNALING")
p4 = plotEnrichment(fgsea_sets[["REACTOME_PI3K_EVENTS_IN_ERBB2_SIGNALING"]],
                    genelist1) + labs(title="REACTOME_PI3K_EVENTS_IN_ERBB2_SIGNALING")

p = (p1 | p4)/(p2|p3)
pdf('fgsea_reactome.pdf',width = 8,height = 6)
p
dev.off()

fwrite(fgseaRes, file="fgseaRes_Reactome.txt", sep="\t", sep2=c("", " ", ""))











