

library(here)
library(tidyverse)
library(vroom)
library(data.table)
library(nflplotR)
library(IOBR)
library(FactoMineR)
library(factoextra)
library(ggridges)
library(ggsci)
library(scales)
library(ggrepel)
library(DESeq2)
library(UpSetR)
library(ComplexUpset)
library(org.Hs.eg.db)
library(GEOmirror)
library(WGCNA)




load('F:/01_Exp_SD.Rdata')

load('F:/01_Exp_Stroke.Rdata')


##1.Differential expression analysis ----



IDs = c('GSE208668','GSE56931','GSE98582')

Exp_list = list(Exp_GSE208668_Array_filter,Exp_GSE56931_Array_filter,Exp_GSE98582_Array_filter)

Pd_list = list(Pd_GSE208668_filter,Pd_GSE56931,Pd_GSE98582)

SD_DEGs_list = list()

for (i in 1:length(IDs)) {
  
  library(limma)
  
  Condition <- Pd_list[[i]]$Group %>% factor(levels = c('Control','SD'))
  
  design = model.matrix(~Condition)   
  
  fit = lmFit(Exp_list[[i]],design)
  
  fit = eBayes(fit)
  
  DEG = topTable(fit,coef=2,number = Inf) %>% 
    as.data.frame() %>% 
    mutate(Change = ifelse(logFC > 0.25 & P.Value < 0.05, 'Up', 
                           ifelse(logFC < -0.25 & P.Value < 0.05, 'Down', 'Not')))
  
  SD_DEGs_list[[i]] = DEG
}


table(SD_DEGs_list[[1]]$Change)
# Down   Not    Up 
# 2948 28072  2189 

table(SD_DEGs_list[[2]]$Change)
# Down   Not    Up 
# 427 20416   261  

table(SD_DEGs_list[[3]]$Change)
# Down  Not   Up 
# 78 6453  226 





p_volcano_GSE208668 <- 
  ggplot(SD_DEGs_list[[1]],aes(x = logFC,y = -log10(P.Value))) +
  geom_point(aes(color = Change),alpha = 0.8,size = 3) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1,
        
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = 'black'),
        axis.title = element_text(color = 'black',size = 10,face = 'bold'),
        title = element_text(color = 'black',size = 12,face = 'bold'),
        panel.grid.major = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_rect(fill = 'white', colour = 'black')) +
  scale_color_manual(name = 'Expression',
                     # color or three types
                     values = c('Up'="#ed6a47",'Not'='grey','Down'="#709fca"),
                     # legend labels
                     label = c('Up'="Up (2189)",'Not'='Not (28072)','Down'='Down (2948)')) +
  # Threshold split line
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 1) +
  geom_vline(xintercept = c(-0.25,0.25),lty = 'dashed',size = 1) + 
  # scale_x_continuous(breaks = c(-4,-2,-1,0,1,2,4)) +
  xlab('Log2FoldChange')+ylab('-Log10(Pvalue)')+
  ggtitle('GSE208668')

ggpreview(p_volcano_GSE208668,width = 5, height = 5)
ggsave(p_volcano_GSE208668,filename= here('',paste0("GSE208668_","volcano.pdf")), width = 5, height = 5)
ggsave(p_volcano_GSE208668,filename= here('',paste0("GSE208668_","volcano.tiff")), width = 5, height = 5, dpi = 300)



#**Enrichment analysis was carried out for DEGs**

Sel_Tra_Enrich_ORA_GSEA <- function(ID = ID, DEGs_dat = DEGs_dat,
                                    minGSSize = minGSSize, maxGSSize = maxGSSize,
                                    fromtype = fromtype, orgdb = orgdb, org = org, pcut = pcut)
{
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  data <- DEGs_dat %>% 
    as.data.frame() %>% 
    filter(Change %in% c('Up','Down')) %>% 
    mutate(gene = rownames(.))
  
  s2e <- bitr(rownames(data), 
              fromType = fromtype,
              toType = "ENTREZID",
              OrgDb = orgdb)
  
  gene_module_entrz <- merge(s2e, data, by.x = fromtype, by.y="gene") %>% 
    arrange(desc(log2FoldChange))
  
  print("Tra finished")
  
  #**ORA**
  
  ego = enrichGO(gene = gene_module_entrz$ENTREZID,
                 OrgDb= orgdb,
                 ont = "ALL",
                 minGSSize = minGSSize,
                 maxGSSize = maxGSSize,
                 pvalueCutoff = pcut,
                 readable = TRUE)
  
  write.csv(ego@result,file = here('Output/1_DEGs_WGCNA_GSEA', paste0('GO_',ID,'_DEGs_.csv')))
  
  print("DEGs_GO finished")
  
  ekegg = enrichKEGG(gene = gene_module_entrz$ENTREZID,
                     pvalueCutoff = pcut,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     organism = org)
  
  write.csv(ekegg@result,file = here('Output/1_DEGs_WGCNA_GSEA', paste0('KEGG_',ID,'_DEGs_.csv')))
  
  print("DEGs_KEGG finished")
  
  save(ego, ekegg, file = here('Data/1_DEGs_WGCNA_GSEA', paste0(ID,'_DEGs_GO_KEGG.Rdata')))
  
  
  
  #**GSEA**
  
  
  gene <- gene_module_entrz$log2FoldChange
  
  names(gene) <- gene_module_entrz$ENTREZID
  
  
  GO <- gseGO(
    gene, 
    ont = "ALL",# "BP"、"MF"和"CC"或"ALL"
    OrgDb = org.Hs.eg.db,#Human annotated genes
    keyType = "ENTREZID",
    pvalueCutoff = pcut,
    pAdjustMethod = "BH"#p-value correction method
  )
  # preparing geneSet collections...
  # GSEA analysis...
  # leading edge analysis...
  # done...
  # There were 14 warnings (use warnings() to see them)
  
  
  sortGO <- GO[order(GO@result$NES, decreasing = T),]
  
  write.csv(sortGO,file = here('Output/1_DEGs_WGCNA_GSEA', paste0('GO_DEGs_GSEA_',ID,'_.csv')))
  
  #KEGG Enrichment analysis
  
  KEGG <- gseKEGG(
    gene,
    organism = "hsa",
    keyType = "kegg",
    exponent = 1,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize ,
    eps = 1e-10,
    pvalueCutoff = pcut,
    pAdjustMethod = "BH",
    verbose = TRUE,
    use_internal_data = FALSE,
    seed = FALSE,
    by = "fgsea")
  # preparing geneSet collections...
  # GSEA analysis...
  # 
  # 
  # leading edge analysis...
  # done...
  
  sortKEGG <- KEGG[order(KEGG@result$NES, decreasing = T),]
  
  write.csv(sortKEGG,file = here('Output/1_DEGs_WGCNA_GSEA', paste0('KEGG_DEGs_GSEA_',ID,'_.csv')))
  
  save(GO, KEGG, sortGO, sortKEGG, file = here('Data/1_DEGs_WGCNA_GSEA', paste0(ID,'_DEGs_GSEA_GO_KEGG.Rdata')))
  
  
}

colnames(SD_DEGs_list[[1]])
# [1] "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"         "Change" 

colnames(SD_DEGs_list[[1]]) = c("log2FoldChange","AveExpr","t","P.Value","adj.P.Val","B","Change")

Sel_Tra_Enrich_ORA_GSEA(ID = 'GSE208668',DEGs_dat = SD_DEGs_list[[1]],
                        minGSSize = 1,maxGSSize = 10000,
                        fromtype = "SYMBOL", orgdb = org.Hs.eg.db, org = 'hsa', pcut = 1)


#Visualize enrichment results

load(here('Data/1_DEGs_WGCNA_GSEA','GSE208668_DEGs_GO_KEGG.Rdata'))

GO_data <- ego@result %>% 
  as.data.frame() %>% 
  filter(pvalue < 0.05 & ONTOLOGY == 'BP') %>% 
  mutate(Description = str_to_title(Description)) %>%
  head(15) %>% 
  arrange(pvalue)


mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), 
  axis.ticks.y = element_blank(), 
  plot.title = element_text(size = 14,
                            hjust = 0.5,
                            face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
                       r = 10,
                       l = 5.5,
                       b = 5.5)
)

library(ggbreak)

p_BP <- ggplot(data = GO_data, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  # scale_x_break(c(45,360))+
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "Number of Gene", y = "", title = "BP of GO enrichment barplot") +
  geom_text(aes(x = 0.03,
                label = Description),
            hjust = 0)+ 
  theme_classic() + mytheme

ggpreview(p_BP,width = 6,height = 6)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_BP_GSE208668_boxplot.pdf'),p_BP,width = 6,height = 6,)
ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_BP_GSE208668_boxplot.tiff'),p_BP,width = 6,height = 6, dpi = 300)




KEGG_data <- ekegg@result %>% 
  as.data.frame() %>% 
  filter(grepl(x = category,'Metabolism') | grepl(x = category,'Organismal Systems')) %>%
  filter(pvalue < 0.05 ) %>% 
  mutate(Description = str_to_title(Description)) %>%
  head(15) %>% 
  arrange(pvalue)

p_KEGG <- ggplot(data = KEGG_data, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  scale_x_break(c(90,340))+
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "Number of Gene", y = "", title = "KEGG enrichment barplot") +
  geom_text(aes(x = 0.03,
                label = Description),
            hjust = 0)+ 
  theme_classic() + mytheme

ggpreview(p_KEGG,width = 6,height = 6)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_KEGG_GSE208668_boxplot.pdf'),p_KEGG,width = 6,height = 6)
ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_KEGG_GSE208668_boxplot.tiff'),p_KEGG,width = 6,height = 6, dpi = 300)













#**Differential gene analysis was performed on the Stroke dataset**

library(limma)

Condition <- Pd_GSE16561$Group %>% factor(levels = c('Control','Stroke'))

design = model.matrix(~Condition)   

fit = lmFit(Exp_GSE16561_Array_filter,design)

fit = eBayes(fit)

DEG = topTable(fit,coef=2,number = Inf) %>% 
  as.data.frame() %>% 
  mutate(Change = ifelse(logFC > 0.25 & P.Value < 0.05, 'Up', 
                         ifelse(logFC < -0.25 & P.Value < 0.05, 'Down', 'Not')))

GSE16561_DEGs = DEG

table(GSE16561_DEGs$Change)
# 
# Down  Not   Up 
# 394 6482  492 


#火山图visualization

p_volcano_GSE16561 <- 
  ggplot(GSE16561_DEGs,aes(x = logFC,y = -log10(P.Value))) +
  geom_point(aes(color = Change),alpha = 0.8,size = 3) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1,
        # 标题居中
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = 'black'),
        axis.title = element_text(color = 'black',size = 10,face = 'bold'),
        title = element_text(color = 'black',size = 12,face = 'bold'),
        panel.grid.major = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_rect(fill = 'white', colour = 'black')) +
  scale_color_manual(name = 'Expression',
                     # color or three types
                     values = c('Up'="#ed6a47",'Not'='grey','Down'="#709fca"),
                     # legend labels
                     label = c('Up'="Up (492)",'Not'='Not (6482)','Down'='Down (394)')) +
  # Threshold split line
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 1) +
  geom_vline(xintercept = c(-0.25,0.25),lty = 'dashed',size = 1) + 
  scale_x_continuous(breaks = c(-2,-1,0,1,2)) +
  xlim(-4,4)+
  xlab('Log2FoldChange')+ylab('-Log10(Pvalue)')+
  ggtitle('GSE16561')

ggpreview(p_volcano_GSE16561,width = 5, height = 5)
ggsave(p_volcano_GSE16561,filename= here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE16561_","volcano.pdf")), width = 5, height = 5)
ggsave(p_volcano_GSE16561,filename= here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE16561_","volcano.tiff")), width = 5, height = 5, dpi = 300)




#**Enrichment analysis was carried out for DEGs**

colnames(GSE16561_DEGs)
# [1] "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"         "Change" 

colnames(GSE16561_DEGs) = c("log2FoldChange","AveExpr","t","P.Value","adj.P.Val","B","Change")

Sel_Tra_Enrich_ORA_GSEA(ID = 'GSE16561',DEGs_dat = GSE16561_DEGs,
                        minGSSize = 1,maxGSSize = 10000,
                        fromtype = "SYMBOL", orgdb = org.Hs.eg.db, org = 'hsa', pcut = 1)


#Visualize enrichment results


load(here('Data/1_DEGs_WGCNA_GSEA','GSE16561_DEGs_GO_KEGG.Rdata'))

GO_data <- ego@result %>% 
  as.data.frame() %>% 
  filter(pvalue < 0.05 & ONTOLOGY == 'BP') %>% 
  #添加命令，使Description首字母大写
  mutate(Description = str_to_title(Description)) %>%
  head(15) %>% 
  arrange(pvalue)

#先自定义主题：
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), 
  axis.ticks.y = element_blank(), 
  plot.title = element_text(size = 14,
                            hjust = 0.5,
                            face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
                       r = 10,
                       l = 5.5,
                       b = 5.5)
)

library(ggbreak)

p_BP <- ggplot(data = GO_data, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  # scale_x_break(c(45,360))+
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "Number of Gene", y = "", title = "BP of GO enrichment barplot") +
  geom_text(aes(x = 0.03,
                label = Description),
            hjust = 0)+ 
  theme_classic() + mytheme

ggpreview(p_BP,width = 6,height = 6)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_BP_GSE16561_boxplot.pdf'),p_BP,width = 6,height = 6,)
ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_BP_GSE16561_boxplot.tiff'),p_BP,width = 6,height = 6, dpi = 300)




KEGG_data <- ekegg@result %>% 
  as.data.frame() %>% 
  filter(grepl(x = category,'Metabolism') | grepl(x = category,'Organismal Systems')) %>%
  filter(pvalue < 0.05 ) %>% 
  mutate(Description = str_to_title(Description)) %>%
  head(15) %>% 
  arrange(pvalue)

p_KEGG <- ggplot(data = KEGG_data, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  scale_x_break(c(20,90))+
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "Number of Gene", y = "", title = "KEGG enrichment barplot") +
  geom_text(aes(x = 0.03, 
                label = Description),
            hjust = 0)+ 
  theme_classic() + mytheme

ggpreview(p_KEGG,width = 6,height = 6)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_KEGG_GSE16561_boxplot.pdf'),p_KEGG,width = 6,height = 6)
ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_KEGG_GSE16561_boxplot.tiff'),p_KEGG,width = 6,height = 6, dpi = 300)







#**Integrated analysis was performed for differentially expressed genes**

#UPSET

Upset_data = fromList(list("UP_DEGs (GSE208668)" = rownames(SD_DEGs_list[[1]][SD_DEGs_list[[1]]$Change =='Up',]),
                           "DOWN_DEGs (GSE208668)" =  rownames(SD_DEGs_list[[1]][SD_DEGs_list[[1]]$Change =='Down',]),
                           "UP_DEGs (GSE56931)" =  rownames(SD_DEGs_list[[2]][SD_DEGs_list[[2]]$Change =='Up',]),
                           "DOWN_DEGs (GSE56931)" =  rownames(SD_DEGs_list[[2]][SD_DEGs_list[[2]]$Change =='Down',]),
                           "UP_DEGs (GSE98582)" = rownames(SD_DEGs_list[[3]][SD_DEGs_list[[3]]$Change =='Up',]),
                           "DOWN_DEGs (GSE98582)" =  rownames(SD_DEGs_list[[3]][SD_DEGs_list[[3]]$Change =='Down',]),
                           "UP_DEGs (GSE16561)" =  rownames(GSE16561_DEGs[GSE16561_DEGs$Change =='Up',]),
                           "DOWN_DEGs (GSE16561)" =  rownames(GSE16561_DEGs[GSE16561_DEGs$Change =='Down',])))

groups = c("UP_DEGs (GSE208668)",
           "DOWN_DEGs (GSE208668)",
           "UP_DEGs (GSE56931)",
           "DOWN_DEGs (GSE56931)",
           "UP_DEGs (GSE98582)",
           "DOWN_DEGs (GSE98582)",
           "UP_DEGs (GSE16561)",
           "DOWN_DEGs (GSE16561)")



p_Upset = ComplexUpset::upset(Upset_data,groups,
                              name = "", 
                              # intersections=list(c("UP_DEGs (Huh7)"),
                              #                    c("DOWN_DEGs (Huh7)"),
                              #                    c("UP_DEGs (Hep3B)"),
                              #                    c("DOWN_DEGs (Hep3B)"),
                              #                    c("UP_DEGs (Huh7)","UP_DEGs (Hep3B)"),
                              #                    c("DOWN_DEGs (Huh7)","DOWN_DEGs (Hep3B)")),                  
                              sort_intersections=FALSE,
                              sort_sets=FALSE,
                              # base_annotations = list(
                              #   "intersection size"=intersection_size(bar_number_threshold = 1,
                              #                                         mode = 'inclusive_intersection')
                              #   + ylab('Intersection size')
                              # ),
                              base_annotations = list("intersection" = intersection_size(bar_number_threshold = 1,
                                                                                         counts = T,
                                                                                         mode = 'inclusive_intersection',
                                                                                         mapping = aes(fill="bars_color"))+
                                                        scale_fill_manual(values = c("bars_color"="skyblue"),guide="none")+
                                                        theme(panel.grid.major = element_blank(),
                                                              panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(),
                                                              text = element_text(size = 14,color = 'black'))+
                                                        ylab('Inclusive Intersection Size')),
                              mode = 'inclusive_intersection',
                              # queries=list(upset_query(intersect=c("UP_DEGs (3_Cluster1)","UP_DEGs (3_Cluster2)"), color="#ed6a47", fill="#ed6a47"),
                              #              upset_query(intersect=c("DOWN_DEGs (3_Cluster1)","DOWN_DEGs (3_Cluster2)"), color="#709fca", fill="#709fca"),
                              #              upset_query(intersect=c("UP_DEGs (3_Cluster1)","DOWN_DEGs (3_Cluster2)"), color='#5F4495', fill='#5F4495'),
                              #              upset_query(intersect=c("DOWN_DEGs (3_Cluster1)","UP_DEGs (3_Cluster2)"), color='#923693', fill='#923693')),
                              # # intersections='all',
                              width_ratio = 0, 
                              height_ratio = 0.55
)+theme(text = element_text(size = 10,color = 'black'))

ggpreview(p_Upset,width = 10,height = 6)



ggsave(p_Upset,width = 10,height = 6,filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_Upset_0.5.pdf'))
ggsave(p_Upset,width = 10,height = 6,filename = here('Output/1_DEGs_WGCNA_GSEA/','DEGs_Upset_0.5.tiff'))







#**Select the test set for VEEN**


Upset_data = fromList(list("UP_DEGs (GSE208668)" = rownames(SD_DEGs_list[[1]][SD_DEGs_list[[1]]$Change =='Up',]),
                           "DOWN_DEGs (GSE208668)" =  rownames(SD_DEGs_list[[1]][SD_DEGs_list[[1]]$Change =='Down',]),
                           "UP_DEGs (GSE16561)" =  rownames(GSE16561_DEGs[GSE16561_DEGs$Change =='Up',]),
                           "DOWN_DEGs (GSE16561)" =  rownames(GSE16561_DEGs[GSE16561_DEGs$Change =='Down',])))

groups = c("UP_DEGs (GSE208668)",
           "DOWN_DEGs (GSE208668)",
           "UP_DEGs (GSE16561)",
           "DOWN_DEGs (GSE16561)")


p_Upset = ComplexUpset::upset(Upset_data,groups,
                              name = "", 
                              intersections=list(c("UP_DEGs (GSE208668)"),
                                                 c("DOWN_DEGs (GSE208668)"),
                                                 c("UP_DEGs (GSE16561)"),
                                                 c("DOWN_DEGs (GSE16561)"),
                                                 c("UP_DEGs (GSE208668)","UP_DEGs (GSE16561)"),
                                                 c("DOWN_DEGs (GSE208668)","DOWN_DEGs (GSE16561)")),
                              sort_intersections=FALSE,
                              sort_sets=FALSE,
                              # base_annotations = list(
                              #   "intersection size"=intersection_size(bar_number_threshold = 1,
                              #                                         mode = 'inclusive_intersection')
                              #   + ylab('Intersection size')
                              # ),
                              base_annotations = list("intersection" = intersection_size(bar_number_threshold = 1,
                                                                                         counts = T,
                                                                                         mode = 'inclusive_intersection',
                                                                                         mapping = aes(fill="bars_color"))+
                                                        scale_fill_manual(values = c("bars_color"="skyblue"),guide="none")+
                                                        theme(panel.grid.major = element_blank(),
                                                              panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(),
                                                              text = element_text(size = 14,color = 'black'))+
                                                        ylab('Inclusive Intersection Size')),
                              mode = 'inclusive_intersection',
                              # queries=list(upset_query(intersect=c("UP_DEGs (3_Cluster1)","UP_DEGs (3_Cluster2)"), color="#ed6a47", fill="#ed6a47"),
                              #              upset_query(intersect=c("DOWN_DEGs (3_Cluster1)","DOWN_DEGs (3_Cluster2)"), color="#709fca", fill="#709fca"),
                              #              upset_query(intersect=c("UP_DEGs (3_Cluster1)","DOWN_DEGs (3_Cluster2)"), color='#5F4495', fill='#5F4495'),
                              #              upset_query(intersect=c("DOWN_DEGs (3_Cluster1)","UP_DEGs (3_Cluster2)"), color='#923693', fill='#923693')),
                              # # intersections='all',
                              width_ratio = 0, 
                              height_ratio = 0.55
)+theme(text = element_text(size = 10,color = 'black'))

ggpreview(p_Upset,width = 6,height = 5)
ggsave(p_Upset,width = 6,height = 5,filename = here('Output/1_DEGs_WGCNA_GSEA/','SD_Stroke_DEGs_Upset.pdf'))
ggsave(p_Upset,width = 6,height = 5,filename = here('Output/1_DEGs_WGCNA_GSEA/','SD_Stroke_DEGs_Upset.tiff'))










##2.SD_WGCNA----

###1. Organize the input data------

library(WGCNA)

enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 



#For RNA-seq data, select log (fpkm+1)
#For chip data, select a uniformized value
data <- as.matrix(Exp_GSE208668_Array_filter)


keep_data <- data[order(apply(data,1,mad), decreasing = T)[1:5000],]

### Create datTraits, which contains grouping, phenotype, and other information
datTraits <- data.frame(row.names = colnames(data),
                        group = Pd_GSE208668_filter$Group)
fix(datTraits)



datExpr0 <- as.data.frame(t(keep_data))



gsg <- goodSamplesGenes(datExpr0,verbose = 3)

gsg$allOK
# [1] TRUE

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste0(names(datExpr0)[!gsg$goodGenes],
                                               collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

gsg <- goodSamplesGenes(datExpr0,verbose = 3)

gsg$allOK
# [1] TRUE



if(T){
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                  signed = FALSE)
  par(mar = c(1,4,3,1),cex=0.8)
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step1_Sample dendrogram and trait.pdf")),width = 8,height = 6)
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )
  ## Plot a line to show the cut
  # abline(h = 23500, col = "red") 
  dev.off()
}


if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 23500, minSize = 10) 
  table(clust)
  keepSamples <- (clust==1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr0) 
}



group_list <- datTraits$group
dat.pca <- PCA(datExpr0, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE, 
                    col.ind = group_list, 
                    axes.linetype=NA,  # remove axeslines
                    mean.point=F
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) 
pca
ggsave(pca,filename= here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step1_Sample PCA analysis.pdf")), width = 8, height = 8)
ggsave(pca,filename= here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step1_Sample PCA analysis.tiff")), width = 8, height = 8, dpi = 300)


##Save data
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file=here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step1_input.Rdata")))



###2. Choose the best threshold power------

# rm(list = ls())  
# 
# load(here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step1_input.Rdata")))

R.sq_cutoff = 0.85  

if(T){
  # Call the network topology analysis function
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  #SFT.R.sq > 0.8 , slope ≈ -1
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step2_power-value.pdf")),width = 16,height = 12)
  
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=R.sq_cutoff ,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  dev.off()
}

sft$powerEstimate  #
# 8

# power = sft$powerEstimate
power = 10

# If the power of the undirected network is less than 15 or the power of the directed network is less than 30, there is no power value
# The scale of the network graph structure R^2 reaches 0.8 and the average connection is below 100, which may be due to the
# Some samples are too different from others. This can be caused by batch effects, sample heterogeneity, or experimental condition pairs
# Too much influence of expression, etc. You can view grouping information and the presence of abnormal samples by plotting sample clustering.
# If this is indeed caused by meaningful biotic changes, you can also use the empirical power values below.



if(is.na(power)){
  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

save(sft, power, file=here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step2_power_value.Rdata")))




###3.A weighted co-expression network was constructed by one-step method to identify gene modules------

# 
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step1_input.Rdata")))
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step2_power_value.Rdata")))

if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr),
    corType = "pearson", 
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 50,    
    mergeCutHeight = 0.15, 
    numericLabels = TRUE, 
    saveTOMs = F,
    verbose = 3
  )
  table(net$colors)
}


if(T){
  # Convert labels to colors for plotting
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step3_genes-modules_ClusterDendrogram.pdf")),width = 16,height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}

save(net, moduleColors, file = here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step3_genes_modules.Rdata")))




###4.Associated gene modules with phenotypic blocks------


rm(list = ls())
load(file = here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step1_input.Rdata")))
load(file = here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step2_power_value.Rdata")))
load(file = here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step3_genes_modules.Rdata")))



if(T){
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0+datTraits$group)
  colnames(design) <- levels(datTraits$group)
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step4_Module-trait-relationship_heatmap.pdf")),
      width = 2*length(colnames(design)), 
      height = 0.6*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3)) #留白：下、左、上、右
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1), 
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = here('Data/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step4_design.Rdata")))
}


if(T){
  mes_group <- merge(MEs,datTraits,by="row.names") 
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  draw_ggboxplot <- function(data,Module="Module",group="group"){
    ggboxplot(data,x=group, y=Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = "jco",
              #add="jitter",
              legend = "") +stat_compare_means()
  }
  # boxplot
  colorNames <- names(MEs)
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step4_Module-trait-relationship_boxplot.pdf")), width = 7.5,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=2)) 
  dev.off()
}






##**5.WGCNA visualization：

### Module relevance display Eigengene-adjacency-heatmap
if(T){
  MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
  MET = orderMEs(MEs)
  # 若添加表型数据
  if(T){
    ## 连续型性状
    # MET = orderMEs(cbind(MEs,datTraits$groupNo))
    design
    SD = as.data.frame(design[,1])
    names(SD) = "SD"
    # Add the weight to existing module eigengenes
    MET = orderMEs(cbind(MEs, SD))
  }
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste0("GSE208668_","step5_module_cor_Eigengene-dendrogram.pdf")),width = 8,height = 10)
  plotEigengeneNetworks(MET, setLabels="", 
                        marDendro = c(0,4,1,4),  
                        marHeatmap = c(5,5,1,2), 
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}





##**6.The gene module of interest was selected for functional enrichment analysis**


##**GO and KEGG analysis is done directly through custom functions**

Sel_Tra_Enrich <- function(ID = ID, exp_dat = exp_dat, module = dynamicColors,
                           minGSSize = minGSSize, maxGSSize = maxGSSize,choose_module = choose_module,
                           fromtype = fromtype, orgdb = orgdb, org = org, pcut = pcut)
{
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gene_module <- data.frame(gene = colnames(exp_dat),
                            module = module)
  
  tmp = bitr(gene_module$gene, fromType = fromtype, toType = 'ENTREZID', OrgDb = orgdb)
  
  print("Tra finished")
  
  gene_module_entrz <- merge(tmp,gene_module, by.x = fromtype, by.y="gene")
  
  choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
  
  ego = enrichGO(gene = choose_gene_module_entrz$ENTREZID,
                 OrgDb= orgdb,
                 ont = "ALL",
                 minGSSize = minGSSize,
                 maxGSSize = maxGSSize,
                 pvalueCutoff = pcut,
                 readable = TRUE)
  write.csv(ego@result,file = here('Output/1_DEGs_WGCNA_GSEA', paste0('GO_',ID,'_',choose_module,'_.csv')))
  print("GO finished")
  ekegg = enrichKEGG(gene = choose_gene_module_entrz$ENTREZID,
                     pvalueCutoff = pcut,
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     organism = org)
  write.csv(ekegg@result,file = here('Output/1_DEGs_WGCNA_GSEA', paste0('KEGG_',ID,'_',choose_module,'_.csv')))
  print("KEGG finished")
  save(ego, ekegg, file = here('Data/1_DEGs_WGCNA_GSEA', paste0(ID,'_',choose_module,'_GO_KEGG.Rdata')))
}


Sel_Tra_Enrich(ID = 'GSE208668',exp_dat = datExpr,module = moduleColors,
               minGSSize = 1,maxGSSize = 10000,choose_module = "blue",
               fromtype = "SYMBOL", orgdb = org.Hs.eg.db, org = 'hsa', pcut = 1)

Sel_Tra_Enrich(ID = 'GSE208668',exp_dat = datExpr,module = moduleColors,
               minGSSize = 1,maxGSSize = 10000,choose_module = "turquoise",
               fromtype = "SYMBOL", orgdb = org.Hs.eg.db, org = 'hsa', pcut = 1)




###5.WGCNA Personalized Visualization------

load(file = here('Data/1_DEGs_WGCNA_GSEA/','GSE208668_step1_input.Rdata'))
load(file = here('Data/1_DEGs_WGCNA_GSEA/','GSE208668_step2_power_value.Rdata'))
load(file = here('Data/1_DEGs_WGCNA_GSEA/','GSE208668_step3_genes_modules.Rdata'))
load(file = here('Data/1_DEGs_WGCNA_GSEA/','GSE208668_step4_design.Rdata'))


####Plot the correlation heatmap between shapes and modules-----


library(tidyverse)
library(reshape2)

moduleTraitCor = cor(net$MEs, design, use = "p");
moduleTraitCor
# Control            SD
# ME2 -0.8758344288  0.8758344288
# ME3 -0.5070662624  0.5070662624
# ME5  0.3305765425 -0.3305765425
# ME1  0.9926282991 -0.9926282991
# ME4  0.6013863355 -0.6013863355
# ME0 -0.0007833164  0.0007833164

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(net$MEs))
moduleTraitPvalue
# Control     SD
# ME2 2.500087e-11 2.500087e-11
# ME3 2.598416e-03 2.598416e-03
# ME5 6.024315e-02 6.024315e-02
# ME1 5.548579e-30 5.548579e-30
# ME4 2.142899e-04 2.142899e-04
# ME0 9.965481e-01 9.965481e-01


trait_heatmap_data <- as.data.frame(moduleTraitCor) %>% 
  mutate(module = rownames(.)) %>% 
  gather(key = cluster,
         value = Cor,
         - module)

trait_heatmap_p <- as.data.frame(moduleTraitPvalue) %>% 
  mutate(module = rownames(.)) %>% 
  gather(key = cluster,
         value = pvalue,
         - module) 

trait_heatmap_data <- trait_heatmap_data %>%
  mutate(pvalue = trait_heatmap_p$pvalue) %>% 
  mutate(text = paste0(round(Cor,3),paste0("\n",paste0('(',format(pvalue,scientific = T,digits = 3),')')))) %>% 
  mutate(cluster = factor(cluster, levels = c("Control", "SD")))




##**左侧注释**


left <- ggplot(data = data.frame(ID = colnames(net$MEs),
                                 y = 'A',
                                 cor = colnames(net$MEs)),
               mapping = aes(y,y=ID,fill=cor))+
  geom_tile() + 
  scale_fill_manual(values = c("#FED439FF","#FD7446FF","#8A9197FF",
                               "#D2AF81FF","#709AE1FF","#D5E4A2FF",
                               "#197EC0FF","#1A9993FF","#FD8CC1FF"))+
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        panel.grid.major=element_blank(),legend.position = 'none')




##**主图**


P_trait_heatmap <- ggplot(trait_heatmap_data, aes(cluster, module)) + 
  geom_tile(aes(fill = Cor), colour = "black", size = 0.5)+
  scale_fill_gradient2(low="#709fca",high="#ed6a47", mid="white", midpoint = 0) + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(hjust = 1, size = 10, face = "bold"),
        panel.grid.major=element_blank()) + 
  labs(fill = "Cor")   


library(aplot)
library(ggplotify)

P_trait_heatmap <- P_trait_heatmap %>% 
  insert_left(left,width = 0.15)

P_trait_heatmap <- as.ggplot(P_trait_heatmap)

P_trait_heatmap <- P_trait_heatmap+
  theme(plot.title = element_text(hjust = 0.6,size = 14,face = 'bold'))+
  labs(title = 'Module-traits relationships')

ggpreview(P_trait_heatmap,width = 4,height = 5)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE208668_WGCNA_Module_Traits.pdf'),P_trait_heatmap,width = 3,height = 5)




####Plot the expression heatmap of genes in the module-----

ME1_Exp <- datExpr %>% 
  t() %>%
  as.data.frame() %>% 
  mutate(ID = net$colors) %>% 
  filter(ID == 1) %>% 
  mutate(ID = NULL)

ME1_Exp[1:5,1:5]
# GSM6360934 GSM6360935 GSM6360936 GSM6360937 GSM6360938
# HS.533945   9.735491   10.51230   9.226758   10.24631   10.94908
# LOC645979  12.209798   12.55280  10.411238   12.19771   10.13576
# LOC441377  12.503275   12.83537  10.508833   12.47516   10.35294
# LOC641768  12.281409   12.80612  10.395440   12.26055   10.30402
# HS.544632  11.028485   10.58990  10.134490   10.72261   10.86414

col_metaData <- data.frame(Group = Pd_GSE208668_filter$Group) %>% 
  mutate(row = colnames(ME1_Exp),
         Group = factor(Group, levels = c("Control", "SD"))) %>% 
  arrange(Group) %>% 
  column_to_rownames('row')

ME1_Exp = ME1_Exp[,match(rownames(col_metaData),colnames(ME1_Exp))]

exprcol <- c("#20854EFF","#E18727FF")
names(exprcol) <- c("Control", "SD")

col <- list(Group=exprcol)

library(ggheatmap)

p <- ggheatmap(ME1_Exp,color = colorRampPalette(c( "#709fca","white","#ed6a47"))(100),
               cluster_rows = T,
               cluster_cols = F,scale = 'row',legendName = 'Nor_Exp',
               text_show_rows = NULL,
               text_show_cols = NULL,
               show_cluster_cols = F,show_cluster_rows = F,
               # border = "brack",
               # cluster_num = c(5,5),
               # annotation_rows = row_metaData,
               annotation_cols = col_metaData,
               annotation_color = col)

ggpreview(p,width = 5.5,height = 5)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE208668_WGCNA_heatmap.pdf'),p,width = 5.5,height = 5)








####Identify genes that may affect phenotype-----

load(file=here('Data/1_DEGs_WGCNA_GSEA/',"GSE208668_step2_power_value.Rdata"))

C1 = design[,"SD",drop=F]

GS_C1 = as.data.frame(cor(datExpr, C1, use = "p"))
colnames(GS_C1) = "GS_C1"
head(GS_C1)
# GS_C1
# HS.533945 0.9518899
# LOC645979 0.7564290
# LOC441377 0.7369184
# LOC641768 0.7467704
# HS.544632 0.9829010
# HS.407822 0.9761419
GS.p_C1 = as.data.frame(corPvalueStudent(as.matrix(GS_C1), nrow(datExpr)))
head(GS.p_C1)
# GS_C1
# HS.533945 1.776399e-17
# LOC645979 3.540742e-07
# LOC441377 1.006565e-06
# LOC641768 6.008872e-07
# HS.544632 2.393416e-24
# HS.407822 3.991460e-22

MM = as.data.frame(cor(datExpr, net$MEs, use = "p"))

point_data <- cbind(MM,GS_C1,GS.p_C1)


#分别绘制不同性状
##**M1**turquoise
module = 'turquoise'
moduleGenes = rownames(MM)[moduleColors=='turquoise']
data <- data.frame(X = abs(point_data[moduleGenes,"ME1"]),
                   Y = abs(point_data[moduleGenes,"GS_C1"]))
P_C1 <- ggplot(data,aes(X,Y))+
  geom_point(size=4,alpha=0.3,color="#FD7446FF")+
  # geom_smooth(method = "lm", formula = y~x, color = "black", fill = "black")+ #颜色选自https://colorbrewer2.org/
  theme_bw()+
  stat_cor()+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', size = 1, fill = 'transparent'),
        plot.title = element_text(hjust = 0.5,size = 14,face = 'bold'))+
  labs(x=paste0("Module Membership in ME1 module"),
       y="Gene significance for SD",
       title = paste0("Module membership vs. Gene significance\n"))

ggpreview(P_C1,width = 5,height = 4)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE208668_WGCNA_scatter.pdf'),P_C1,width = 5,height = 4)









##3.Stroke_WGCNA----

###1. Organize the input data------

library(WGCNA)

### 
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 



#For RNA-seq data, select log (fpkm+1)
#For chip data, select a uniformized value
data <- as.matrix(Exp_GSE16561_Array_filter)

keep_data <- data[order(apply(data,1,mad), decreasing = T)[1:5000],]

### Create datTraits, which contains grouping, phenotype, and other information
datTraits <- data.frame(row.names = colnames(data),
                        group = Pd_GSE16561$Group)
fix(datTraits)



datExpr0 <- as.data.frame(t(keep_data))


gsg <- goodSamplesGenes(datExpr0,verbose = 3)

gsg$allOK
# [1] TRUE

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

gsg <- goodSamplesGenes(datExpr0,verbose = 3)

gsg$allOK
# [1] TRUE



if(T){
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                  signed = FALSE)
  par(mar = c(1,4,3,1),cex=0.8)
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step1_Sample dendrogram and trait.pdf")),width = 8,height = 6)
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )
  dev.off()
}


if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 23500, minSize = 10) 
  table(clust)
  keepSamples <- (clust==1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr0) 
}



group_list <- datTraits$group
dat.pca <- PCA(datExpr0, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE, 
                    col.ind = group_list, 
                    axes.linetype=NA,  # remove axeslines
                    mean.point=F
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) 
pca
ggsave(pca,filename= here('Output/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step1_Sample PCA analysis.pdf")), width = 8, height = 8)
ggsave(pca,filename= here('Output/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step1_Sample PCA analysis.tiff")), width = 8, height = 8, dpi = 300)


##Save data
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file=here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step1_input.Rdata")))



###2. Choose the best threshold power------

# rm(list = ls())  
# 
# load(here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step1_input.Rdata")))

R.sq_cutoff = 0.9  #设置R^2 cut-off，默认为0.85

if(T){
  # Call the network topology analysis function
  #设置power参数选择范围
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  #SFT.R.sq > 0.8 , slope ≈ -1
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step2_power-value.pdf")),width = 16,height = 12)
  # Plot the results: 寻找拐点，确认最终power取值
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=R.sq_cutoff ,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  dev.off()
}

sft$powerEstimate  #查看估计的最佳power
# 12

# power = sft$powerEstimate
power = 10




if(is.na(power)){
  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

save(sft, power, file=here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step2_power_value.Rdata")))




###3.A weighted co-expression network was constructed by one-step method to identify gene modules------
# 
# 
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step1_input.Rdata")))
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step2_power_value.Rdata")))

if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr),
    corType = "pearson", 
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 50,    
    mergeCutHeight = 0.15, 
    numericLabels = TRUE, 
    saveTOMs = F,
    verbose = 3
  )
  table(net$colors) 
}


if(T){
  # Convert labels to colors for plotting
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step3_genes-modules_ClusterDendrogram.pdf")),width = 16,height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}

save(net, moduleColors, file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step3_genes_modules.Rdata")))




###4.Associated gene modules with phenotypic blocks------

# 
# rm(list = ls())  
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step1_input.Rdata")))
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step2_power_value.Rdata")))
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step3_genes_modules.Rdata")))
# 



if(T){
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0+datTraits$group)
  colnames(design) <- levels(datTraits$group)
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step4_Module-trait-relationship_heatmap.pdf")),
      width = 2*length(colnames(design)), 
      height = 0.6*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3)) 
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1), 
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step4_design.Rdata")))
}


 
if(T){
  mes_group <- merge(MEs,datTraits,by="row.names") 
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  draw_ggboxplot <- function(data,Module="Module",group="group"){
    ggboxplot(data,x=group, y=Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = "jco",
              #add="jitter",
              legend = "") +stat_compare_means()
  }
  # 批量画boxplot
  colorNames <- names(MEs)
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step4_Module-trait-relationship_boxplot.pdf")), width = 7.5,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=2)) 
  dev.off()
}





##**5.WGCNA visualization：TOMplot  Eigengene-adjacency-heatmap**


# rm(list = ls())  
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_",'step1_input.Rdata')))
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step2_power_value.Rdata")))
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step3_genes_modules.Rdata")))
# load(file = here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step4_design.Rdata")))

if(T){
  TOM=TOMsimilarityFromExpr(datExpr,power=power)
  dissTOM=1-TOM
  ## draw all genes 
  if(T){
    geneTree = net$dendrograms[[1]]
    plotTOM = dissTOM^7
    diag(plotTOM)=NA
    png(here('Output/1_DEGs_WGCNA_GSEA',paste("GSE20866","step5_TOMplot_Network-heatmap.pdf")),width = 800, height=600) 
    TOMplot(plotTOM,geneTree,moduleColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot")
    dev.off()
  }
  ### draw selected genes to save time...just for test...
  if(F){
    nSelect =0.1*nGenes
    set.seed(123)
    select=sample(nGenes,size = nSelect)
    selectTOM = dissTOM[select,select]
    selectTree = hclust(as.dist(selectTOM),method = "average")
    selectColors = moduleColors[select]
    plotDiss=selectTOM^7
    diag(plotDiss)=NA
    pdf(here('Output/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step5_select_TOMplot_Network-heatmap.pdf")),width=8, height=6)
    TOMplot(plotDiss,selectTree,selectColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot of selected gene")
    dev.off()
  }
}


### Module relevance display Eigengene-adjacency-heatmap
if(T){
  MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
  MET = orderMEs(MEs)
  # 若添加表型数据
  if(T){
    ## 连续型性状
    # MET = orderMEs(cbind(MEs,datTraits$groupNo))
    ## 非连续型性状，需将是否属于这个表型进行0,1数值化，已存于design中
    design
    Insomnia = as.data.frame(design[,1])
    names(Insomnia) = "Insomnia"
    # Add the weight to existing module eigengenes
    MET = orderMEs(cbind(MEs, Insomnia))
  }
  pdf(here('Output/1_DEGs_WGCNA_GSEA',paste("GSE16561_","step5_module_cor_Eigengene-dendrogram.pdf")),width = 8,height = 10)
  plotEigengeneNetworks(MET, setLabels="", 
                        marDendro = c(0,4,1,4),  # 留白：下右上左
                        marHeatmap = c(5,5,1,2), # 留白：下右上左
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}





##**6.The gene module of interest was selected for functional enrichment analysis**


##**GO and KEGG analysis is done directly through custom functions**

Sel_Tra_Enrich(ID = 'GSE16561',exp_dat = datExpr,module = moduleColors,
               minGSSize = 1,maxGSSize = 10000,choose_module = "blue",
               fromtype = "SYMBOL", orgdb = org.Hs.eg.db, org = 'hsa', pcut = 1)




###5.WGCNA Personalized Visualization------


load(file = here('Data/1_DEGs_WGCNA_GSEA/','GSE16561_step1_input.Rdata'))
load(file = here('Data/1_DEGs_WGCNA_GSEA/','GSE16561_step2_power_value.Rdata'))
load(file = here('Data/1_DEGs_WGCNA_GSEA/','GSE16561_step3_genes_modules.Rdata'))
load(file = here('Data/1_DEGs_WGCNA_GSEA/','GSE16561_step4_design.Rdata'))



####Plot the correlation heatmap between shapes and modules-----


library(tidyverse)
library(reshape2)

moduleTraitCor = cor(net$MEs, design, use = "p");
moduleTraitCor
# Control       Stroke
# ME1  -0.024896590  0.024896590
# ME8  -0.473994001  0.473994001
# ME3   0.159488271 -0.159488271
# ME7   0.081768394 -0.081768394
# ME2   0.561334013 -0.561334013
# ME4   0.482770405 -0.482770405
# ME5   0.555937674 -0.555937674
# ME6  -0.560548836  0.560548836
# ME9  -0.112469776  0.112469776
# ME10  0.046757395 -0.046757395
# ME0   0.004039044 -0.004039044

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(net$MEs))
moduleTraitPvalue
# Control       Stroke
# ME1  8.464238e-01 8.464238e-01
# ME8  8.705871e-05 8.705871e-05
# ME3  2.118263e-01 2.118263e-01
# ME7  5.240674e-01 5.240674e-01
# ME2  1.702472e-06 1.702472e-06
# ME4  6.144685e-05 6.144685e-05
# ME5  2.243884e-06 2.243884e-06
# ME6  1.772801e-06 1.772801e-06
# ME9  3.801532e-01 3.801532e-01
# ME10 7.159381e-01 7.159381e-01
# ME0  9.749369e-01 9.749369e-01


trait_heatmap_data <- as.data.frame(moduleTraitCor) %>% 
  mutate(module = rownames(.)) %>% 
  gather(key = cluster,
         value = Cor,
         - module)

trait_heatmap_p <- as.data.frame(moduleTraitPvalue) %>% 
  mutate(module = rownames(.)) %>% 
  gather(key = cluster,
         value = pvalue,
         - module) 

trait_heatmap_data <- trait_heatmap_data %>%
  mutate(pvalue = trait_heatmap_p$pvalue) %>% 
  mutate(text = paste0(round(Cor,3),paste0("\n",paste0('(',format(pvalue,scientific = T,digits = 3),')')))) %>% 
  mutate(cluster = factor(cluster, levels = c("Control", "Stroke")),
         module = factor(module, levels = c("ME0","ME1","ME2","ME3","ME4","ME5","ME6","ME7","ME8","ME9","ME10","ME11")))


##**左侧注释**

left <- ggplot(data = data.frame(ID = colnames(net$MEs),
                                 y = 'A',
                                 cor = colnames(net$MEs)),
               mapping = aes(y,y=ID,fill=cor))+
  geom_tile() + 
  scale_fill_manual(values = c("#FED439FF","#FD7446FF","#8A9197FF",
                               "#D2AF81FF","#709AE1FF","#D5E4A2FF",
                               "#197EC0FF","#1A9993FF","#FD8CC1FF",
                               "#3182bd","#756bb1"))+
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        panel.grid.major=element_blank(),legend.position = 'none')




##**主图**


P_trait_heatmap <- ggplot(trait_heatmap_data, aes(cluster, module)) + 
  geom_tile(aes(fill = Cor), colour = "black", size = 0.5)+
  scale_fill_gradient2(low="#709fca",high="#ed6a47", mid="white", midpoint = 0) + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(hjust = 1, size = 10, face = "bold"),
        panel.grid.major=element_blank()) + 
  labs(fill = "Cor")   


library(aplot)
library(ggplotify)

P_trait_heatmap <- P_trait_heatmap %>% 
  insert_left(left,width = 0.15)

P_trait_heatmap <- as.ggplot(P_trait_heatmap)

P_trait_heatmap <- P_trait_heatmap+
  theme(plot.title = element_text(hjust = 0.6,size = 14,face = 'bold'))+
  labs(title = 'Module-traits relationships')

ggpreview(P_trait_heatmap,width = 3,height = 8)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE16561_WGCNA_Module_Traits.pdf'),P_trait_heatmap,width = 3,height = 8)




####Plot the expression heatmap of genes in the module-----

ME2_Exp <- datExpr %>% 
  t() %>%
  as.data.frame() %>% 
  mutate(ID = net$colors) %>% 
  filter(ID == 2) %>% 
  mutate(ID = NULL)

ME2_Exp[1:5,1:5]
# GSM416528 GSM416529 GSM416530 GSM416531 GSM416532
# GRAP     0.9542662 0.4868613 0.7271056 0.8810299 0.6433352
# PASK     0.7016681 0.8023100 0.8375435 0.8394344 1.0831027
# SPOCK2   1.4691020 0.6367597 0.8257047 0.8985197 1.1281400
# PCED1B   1.3246155 0.6213100 0.8359412 0.8923369 1.0415122
# SLC25A42 1.4794002 0.6797515 0.8027027 0.6929477 0.8548015

col_metaData <- data.frame(Group = Pd_GSE16561$Group) %>% 
  mutate(row = colnames(ME2_Exp)) %>% 
  column_to_rownames('row')

exprcol <- c("#008EA0FF","#C71000FF")
names(exprcol) <- c("Control", "Stroke")

col <- list(Group=exprcol)

library(ggheatmap)

p <- ggheatmap(ME2_Exp,color = colorRampPalette(c( "#709fca","white","#ed6a47"))(100),
               cluster_rows = F,scale = 'row',
               cluster_cols = F,legendName = 'Nor_Exp',
               text_show_rows = NULL,
               text_show_cols = NULL,
               # border = "brack",
               # cluster_num = c(5,5),
               # annotation_rows = row_metaData,
               annotation_cols = col_metaData,
               annotation_color = col)

ggpreview(p,width = 5.5,height = 5)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE16561_WGCNA_heatmap.pdf'),p,width = 5.5,height = 5)








####Identify genes that may affect phenotype-----

# load(file=here('Data/1_DEGs_WGCNA_GSEA',paste("GSE16561_step2_power_value.Rdata")))

C2 = design[,"Stroke",drop=F]
# (1) Gene significance，GS：即比较样本某个基因与对应表型的相关性
GS_C2 = as.data.frame(cor(datExpr, C2, use = "p"))
colnames(GS_C2) = "GS_C2"
head(GS_C2)
# GS_C2
# SLC4A1  0.0729580602
# KRT1   -0.0302025039
# FTH1P3  0.3159934646
# PCDHB9 -0.0004499526
# IFI27   0.1860427804
# ARG1    0.3335043762
GS.p_C2 = as.data.frame(corPvalueStudent(as.matrix(GS_C2), nrow(datExpr)))
head(GS.p_C2)
# GS_C2
# SLC4A1 0.569866163
# KRT1   0.814225782
# FTH1P3 0.011639096
# PCDHB9 0.997207513
# IFI27  0.144323684
# ARG1   0.007561642

MM = as.data.frame(cor(datExpr, net$MEs, use = "p"))

point_data <- cbind(MM,GS_C2,GS.p_C2)



##**M2**red
module = 'blue'
moduleGenes = rownames(MM)[moduleColors=='blue']
data <- data.frame(X = abs(point_data[moduleGenes,"ME2"]),
                   Y = abs(point_data[moduleGenes,"GS_C2"]))
P_C2 <- ggplot(data,aes(X,Y))+
  geom_point(size=4,alpha=0.3,color="#D2AF81FF")+
  # geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+ #颜色选自https://colorbrewer2.org/
  theme_bw()+
  stat_cor()+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', size = 1, fill = 'transparent'),
        plot.title = element_text(hjust = 0.5,size = 14,face = 'bold'))+
  labs(x=paste("Module Membership in ME2 module"),
       y="Gene significance for Stroke",
       title = paste("Module membership vs. Gene significance\n"))

ggpreview(P_C2,width = 5,height = 4)

ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE16561_WGCNA_scatter.pdf'),P_C2,width = 5,height = 4)
















##4.Identify co-enrichment pathways, differential genes, etc----




#**Common enrichment pathways**


load(here('Data/1_DEGs_WGCNA_GSEA/','GSE208668_turquoise_GO_KEGG.Rdata'))

GSE208668_turquoise_GO <- ego@result %>% 
  as.data.frame() %>%
  filter(pvalue < 0.05) %>%
  filter(ONTOLOGY == 'BP') %>% 
  arrange(pvalue) %>% 
  mutate(Group = 'GSE208668')

GSE208668_turquoise_KEGG <- ekegg@result %>% 
  as.data.frame() %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) %>% 
  mutate(Group = 'GSE208668')


load(here('Data/1_DEGs_WGCNA_GSEA/','GSE16561_blue_GO_KEGG.Rdata'))

GSE16561_blue_GO <- ego@result %>% 
  as.data.frame() %>%
  filter(pvalue < 0.05) %>%
  filter(ONTOLOGY == 'BP') %>% 
  arrange(pvalue) %>% 
  mutate(Group = 'GSE16561')

GSE16561_blue_KEGG <- ekegg@result %>% 
  as.data.frame() %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) %>% 
  mutate(Group = 'GSE16561')

table(GSE208668_turquoise_GO$Description %in% GSE16561_blue_GO$Description)
# FALSE  TRUE 
# 1067    60

table(GSE208668_turquoise_KEGG$Description %in% GSE16561_blue_KEGG$Description)
# FALSE  TRUE 
# 69     4 




GO_data = rbind(GSE208668_turquoise_GO,GSE16561_blue_GO) %>% 
  filter(Description %in% GSE208668_turquoise_GO$Description[GSE208668_turquoise_GO$Description %in% GSE16561_blue_GO$Description]) %>%
  filter(ONTOLOGY == 'BP') %>%
  mutate(Group = factor(Group,levels = c('GSE208668','GSE16561')),
         Label = ifelse(pvalue>0.05,'',ifelse(pvalue<=0.05&pvalue>0.01,'*',
                                              ifelse(pvalue<=0.01&pvalue>0.001,'**','***')))) %>% 
  arrange(pvalue)

GO_data = GO_data %>% 
  filter(Description %in% head(GO_data$Description,15)) %>% 
  
  mutate(Description = stringr::str_to_title(Description))

GO_heatmap <-
  ggplot(GO_data,aes(x=Group,y=Description,fill=zScore))+ 
  scale_y_discrete(position = 'right')+
  scale_fill_distiller(palette = "Blues", direction = 1) +
  # scale_fill_gradient2(low="white", high="#ed6a47")+
  # scale_fill_continuous(guide = guide_legend())+
  geom_tile(width = 1,
            height = 1)+
  theme_minimal()+
  geom_text(aes(label=Label),col ="black",size = 8) +
  theme(axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5,size = 12),
        axis.text.y =element_text(size = 12,face = 'bold'),
        legend.position = 'bottom',
        panel.grid.major=element_blank())+
  geom_vline(xintercept = c(1.5,2.5),size=.8,color = 'white')+
  geom_hline(yintercept = c(seq(from = 1.5,to = 19.5,by = 1)),size=.8,color = 'white')+
  xlab(NULL) + ylab(NULL)

ggpreview(GO_heatmap,width = 6,height = 7.5)
ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE16561_WGCNA_GO_heatmap.pdf'),GO_heatmap,width = 6,height = 7.5)
ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE16561_WGCNA_GO_heatmap.tiff'),GO_heatmap,width = 6,height = 7.5,dpi = 300)




KEGG_data = rbind(GSE208668_turquoise_KEGG,GSE16561_blue_KEGG) %>%
  filter(Description %in% GSE208668_turquoise_KEGG$Description[GSE208668_turquoise_KEGG$Description %in% GSE16561_blue_KEGG$Description]) %>%
  arrange(pvalue,Description) %>% 
  mutate(Group = factor(Group,levels = c('GSE208668','GSE16561')),
         Label = ifelse(pvalue>0.05,'',ifelse(pvalue<=0.05&pvalue>0.01,'*',
                                              ifelse(pvalue<=0.01&pvalue>0.001,'**','***')))) %>% 
  arrange(pvalue) %>% 
  
  mutate(Description = stringr::str_to_title(Description))

KEGG_heatmap <-
  ggplot(KEGG_data,aes(x=Group,y=Description,fill=zScore))+ 
  scale_y_discrete(position = 'right')+
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  # scale_fill_continuous(guide = guide_legend())+
  geom_tile(width = 1,
            height = 1)+
  theme_minimal()+
  geom_text(aes(label=Label),col ="black",size = 8) +
  theme(axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5,size = 12),
        axis.text.y =element_text(size = 12,face = 'bold'),
        legend.position = 'bottom',
        panel.grid.major=element_blank())+
  geom_vline(xintercept = c(1.5,2.5),size=.8,color = 'white')+
  geom_hline(yintercept = c(seq(from = 1.5,to = 19.5,by = 1)),size=.8,color = 'white')+
  xlab(NULL) + ylab(NULL)

ggpreview(KEGG_heatmap,width = 4,height = 4.5)
ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE16561_WGCNA_KEGG_heatmap.pdf'),KEGG_heatmap,width = 4,height = 4.5)
ggsave(filename = here('Output/1_DEGs_WGCNA_GSEA/','GSE16561_WGCNA_KEGG_heatmap.tiff'),KEGG_heatmap,width = 4,height = 4.5,dpi = 300)








#**common genes**


load(here('Data/1_DEGs_WGCNA_GSEA/','GSE208668_step3_genes_modules.Rdata'))

GSE208668_turquoise_genes <- moduleColors %>% 
  as.data.frame() %>% 
  mutate(Genes = names(net$colors)) %>% 
  dplyr::rename('Colors' = ".") %>% 
  dplyr::select(Genes,Colors) %>%
  filter(Colors == 'turquoise')


load(here('Data/1_DEGs_WGCNA_GSEA/','GSE16561_step3_genes_modules.Rdata'))

GSE16561_blue_genes <- moduleColors %>% 
  as.data.frame() %>% 
  mutate(Genes = names(net$colors)) %>% 
  dplyr::rename('Colors' = ".") %>% 
  dplyr::select(Genes,Colors) %>%
  filter(Colors == 'blue')

table(GSE208668_turquoise_genes$Genes %in% GSE16561_blue_genes$Genes)
# FALSE  TRUE 
# 3097    61 


#The results of WGCNA and DEGs were combined for analysis

GSE208668_DEGs = SD_DEGs_list[[1]]

library(ggvenn)
xx <- list(GSE208668_ME1 = GSE208668_turquoise_genes$Genes,
           GSE208668_DEGs = rownames(GSE208668_DEGs[!GSE208668_DEGs$Change=='Not',]),
           GSE16561_ME2 = GSE16561_blue_genes$Genes,
           GSE16561_DEGs = rownames(GSE16561_DEGs[!GSE16561_DEGs$Change=='Not',]))

p1 <- ggvenn(xx,show_percentage = T,show_elements = F,label_sep = ",",
             digits = 1,stroke_color = "white",
             text_color = "black",text_size = 8,set_name_size = 8,
             fill_color = c("#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"),
             set_name_color = c("#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))

ggpreview(p1,width = 8,height = 8)

ggsave(p1,width = 8,height = 8,filename = here('Output/1_DEGs_WGCNA_GSEA/','WGCNA_DEGs_VENN.pdf'))
ggsave(p1,width = 8,height = 8,filename = here('Output/1_DEGs_WGCNA_GSEA/','WGCNA_DEGs_VENN.tiff'),dpi = 300)



library(clusterProfiler)
library(org.Hs.eg.db)

gene_module <- intersect(GSE208668_turquoise_genes$Genes,GSE16561_blue_genes$Genes)

tmp = bitr(gene_module, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

ego = enrichGO(gene = tmp$ENTREZID,
               OrgDb= org.Hs.eg.db,
               ont = "ALL",
               minGSSize = 1,
               maxGSSize = 10000,
               pvalueCutoff = 1,
               readable = TRUE)

ekegg = enrichKEGG(gene = tmp$ENTREZID,
                   minGSSize = 1,
                   maxGSSize = 10000,
                   pvalueCutoff = 1,
                   organism = 'hsa')

GO_dat = ego@result %>% 
  as.data.frame() %>%
  filter(pvalue < 0.05) %>%
  filter(ONTOLOGY == 'BP') %>% 
  arrange(pvalue) %>% 
  mutate(Description = stringr::str_to_title(Description)) %>% 
  head(15)

KEGG_dat = ekegg@result %>% 
  as.data.frame() %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) %>% 
  mutate(Description = stringr::str_to_title(Description))




mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11), 
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
                 legend.title = element_text(size = 13), 
                 legend.text = element_text(size = 11)) 

p <- ggplot(data = GO_dat, 
            aes(x = Count, y = Description)) + 
  geom_point(aes(size = Count, color = -log10(pvalue))) + 
  scale_color_distiller(palette = "BrBG",direction = -1) +
  labs(x = "Gene Number", 
       y = "",
       title = "Dotplot of Enriched Biological Process",
       size = "gene number") + 
  theme_bw() +
  mytheme

ggpreview(p,width = 9,height = 6)
ggsave(p,width = 9,height = 6,filename = here('Output/1_DEGs_WGCNA_GSEA/','Int_Genes_GO.pdf'))
ggsave(p,width = 9,height = 6,filename = here('Output/1_DEGs_WGCNA_GSEA/','Int_Genes_GO.tiff'),dpi = 300)



p <- ggplot(data = KEGG_dat, 
            aes(x = Count, y = Description)) + 
  geom_point(aes(size = Count, color = -log10(pvalue))) + 
  scale_color_distiller(palette = "Spectral",direction = -1) +
  labs(x = "Gene Number", 
       y = "",
       title = "Dotplot of Enriched KEGG Pathways",
       size = "gene number") + 
  theme_bw() +
  mytheme

ggpreview(p,width = 6,height = 5)
ggsave(p,width = 6,height = 5,filename = here('Output/1_DEGs_WGCNA_GSEA/','Int_Genes_KEGG.pdf'))
ggsave(p,width = 6,height = 5,filename = here('Output/1_DEGs_WGCNA_GSEA/','Int_Genes_KEGG.tiff'),dpi = 300)
