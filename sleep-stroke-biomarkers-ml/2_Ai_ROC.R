


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
library(shapviz)



##1.Identify characteristic genes based on multiple machine learning algorithms----

#载入WGCNA结果
#load(GSE208668_step3_genes_modules.Rdata')

GSE208668_turquoise_genes <- moduleColors %>% 
  as.data.frame() %>% 
  mutate(Genes = names(net$colors)) %>% 
  dplyr::rename('Colors' = ".") %>% 
  dplyr::select(Genes,Colors) %>%
  filter(Colors == 'turquoise')


#load('GSE16561_step3_genes_modules.Rdata')

GSE16561_blue_genes <- moduleColors %>% 
  as.data.frame() %>% 
  mutate(Genes = names(net$colors)) %>% 
  dplyr::rename('Colors' = ".") %>% 
  dplyr::select(Genes,Colors) %>%
  filter(Colors == 'blue')

DEGs_GSE208668 = rownames(SD_DEGs_list[[1]][!SD_DEGs_list[[1]]$Change =='Not',])

DEGs_GSE16561 =  rownames(GSE16561_DEGs[!GSE16561_DEGs$Change =='Not',])

Genes = intersect(intersect(DEGs_GSE208668,DEGs_GSE16561),
                  intersect(GSE208668_turquoise_genes$Genes,GSE16561_blue_genes$Genes))




###Characteristic gene identification was performed on SD data----

#load('01_Exp_SD.Rdata')

table(colnames(Exp_GSE208668_Array_filter)==rownames(Pd_GSE208668_filter))
# TRUE 
# 33 


library(readr)
library(VIM)
library(caret)
library(rpart)
library(rpart.plot)
library(Metrics)
library(stringr)
library(rpart)
library(tibble)
library(bitops)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(tidyverse)
library(limma)
library(pheatmap)
library(visNetwork)
library(ggpol)
library(ggplot2)
library(sparkline)
library(randomForest)
library(venn)
library(sparkline)
library(dplyr)
library(tidyverse)
library(caret)
library(DALEX)
library(gbm)
library(caret)
library(glmnet) 
library(xgboost)
library(DALEX)
library(gbm)
library(VennDiagram)
library(limma)  
library(neuralnet)
library(NeuralNetTools)



##**GSE208668数据**


data = Exp_GSE208668_Array_filter %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'Gene') %>%
  filter(Gene %in% Genes) %>% 
  column_to_rownames(var = 'Gene')

data = t(data)

group = as.character(Pd_GSE208668_filter$Group)#样本分型



#####RandomForest#####
set.seed(123456)
rf  =  randomForest(as.factor(group)~., data = data, ntree = 500)
optionTrees  =  which.min(rf$err.rate[,1])
rf2 = randomForest(as.factor(group)~., data = data, ntree = optionTrees)
importance = importance(x = rf2)
rfGenes = importance[order(importance[,"MeanDecreaseGini"], decreasing  =  TRUE),]
RF_GSE208668_data = data.frame(Genes = names(rfGenes), 
                               Importance = rfGenes,
                               Method = 'RandomForest')



#####Lasso#####
set.seed(123)
x = data
y = group
fit = glmnet(x, y, family  =  "binomial", alpha = 1)
cvfit = cv.glmnet(x, y, family = "binomial", alpha = 1,type.measure = 'deviance',nfolds  =  10)
coef = coef(fit, s  =  cvfit$lambda.min)
index = which(coef != 0)
lassoGene = row.names(coef)[index]
lassoGenes = lassoGene[-1]
dat = as.matrix(coef)
Lasso_GSE208668_data = data.frame(Genes = rownames(dat)[-1],
                                  Importance = abs(dat[,1])[-1],
                                  Method = 'Lasso')

#####XGboost#####
TrainControl <- trainControl( method  =  "repeatedcv", number  =  10, repeats  =  4)
model<- train(x = data,y = as.factor(group),  method  =  "xgbTree", trControl  =  TrainControl,verbose  =  FALSE)
importance <- varImp(model)
XGBimportant <- as.matrix(importance$importance) 
XGBGenes = XGBimportant[order(XGBimportant[,"Overall"], decreasing  =  TRUE),]
# XGBGenes = names(XGBGenes[XGBGenes>importantvlue])
XGB_GSE208668_data = data.frame(Genes = names(XGBGenes),
                                Importance = XGBGenes,
                                Method = 'Xgboost')



#####Summarize multiple results of calculations#####
AI_GSE208668_data = rbind(RF_GSE208668_data,
                          Lasso_GSE208668_data,
                          XGB_GSE208668_data) %>% 
  mutate(Method = factor(Method, levels = c('RandomForest','Lasso','Xgboost')),
         Genes = factor(Genes))


library(ggsci)

p <- ggplot(data=AI_GSE208668_data, 
            aes(x = Genes, y = Importance, fill = Genes)) +
  geom_bar(stat="identity") +
  #利用ggsci包中的调色板
  scale_fill_manual(values = df$color) +
  labs(x = "",title = 'GSE208668') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1)) +
  facet_wrap(~Method, scales = "free_y",ncol = 1)

ggpreview(p,width = 5,height = 6)
ggsave('AI_GSE208668_data.pdf',p, width = 5,height = 6)
ggsave('AI_GSE208668_data.tiff',p, width = 5,height = 6,dpi = 300)







##**GSE16561数据**


data = Exp_GSE16561_Array_filter %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'Gene') %>%
  filter(Gene %in% Genes) %>% 
  column_to_rownames(var = 'Gene')

data = t(data)

group = as.character(Pd_GSE16561$Group)#样本分型




#####RandomForest#####
set.seed(123456)
rf  =  randomForest(as.factor(group)~., data = data, ntree = 500)
optionTrees  =  which.min(rf$err.rate[,1])
rf2 = randomForest(as.factor(group)~., data = data, ntree = optionTrees)
importance = importance(x = rf2)
rfGenes = importance[order(importance[,"MeanDecreaseGini"], decreasing  =  TRUE),]
RF_GSE16561_data = data.frame(Genes = names(rfGenes), 
                              Importance = rfGenes,
                              Method = 'RandomForest')



#####Lasso#####
set.seed(123)
x = data
y = group
fit = glmnet(x, y, family  =  "binomial", alpha = 1)
cvfit = cv.glmnet(x, y, family = "binomial", alpha = 1,type.measure = 'deviance',nfolds  =  10)
coef = coef(fit, s  =  cvfit$lambda.min)
index = which(coef != 0)
lassoGene = row.names(coef)[index]
lassoGenes = lassoGene[-1]
dat = as.matrix(coef)
Lasso_GSE16561_data = data.frame(Genes = rownames(dat)[-1],
                                 Importance = abs(dat[,1])[-1],
                                 Method = 'Lasso')

#####XGboost#####
TrainControl <- trainControl( method  =  "repeatedcv", number  =  10, repeats  =  4)
model<- train(x = data,y = as.factor(group),  method  =  "xgbTree", trControl  =  TrainControl,verbose  =  FALSE)
importance <- varImp(model)
XGBimportant <- as.matrix(importance$importance) 
XGBGenes = XGBimportant[order(XGBimportant[,"Overall"], decreasing  =  TRUE),]
XGB_GSE16561_data = data.frame(Genes = names(XGBGenes),
                               Importance = XGBGenes,
                               Method = 'Xgboost')



#####Summarize multiple results of calculations总#####
AI_GSE16561_data = rbind(RF_GSE16561_data,
                         Lasso_GSE16561_data,
                         XGB_GSE16561_data) %>% 
  mutate(Method = factor(Method, levels = c('RandomForest','Lasso','Xgboost')),
         Genes = factor(Genes))

#构建离散颜色变量

library(tidyverse)
library(scales)
library(ggsci)
library(magrittr)
library(wesanderson)
library(MetBrewer)

sci_palettes <- list(aaas=pal_aaas()(10),
                     npg=pal_npg()(10),
                     nejm=pal_nejm()(8),
                     lancet=pal_lancet()(9),
                     jama=pal_jama()(7),
                     jco=pal_jco()(10),
                     d3=pal_d3()(10),
                     locuszoom=pal_locuszoom()(7),
                     uchicago=pal_uchicago()(9),
                     startek=pal_startrek()(7),
                     tron=pal_tron()(7),
                     futurama=pal_futurama()(12),
                     simpsons=pal_simpsons()(16),
                     cosmic=pal_cosmic("hallmarks_light")(10),
                     rickandmorty=pal_rickandmorty("schwifty")(12),
                     flatui=pal_flatui("default")(10),
                     frontiers=pal_frontiers("default")(10))


df <- unname(unlist(sci_palettes)) %>% as.data.frame() %>%
  bind_rows(unlist(lapply(names(wes_palettes), wes_palette)) %>%
              as.data.frame()) %>%
  bind_rows(unlist(lapply(names(MetBrewer::MetPalettes), met.brewer)) %>%
              as.data.frame()) %>%
  distinct() %>%
  mutate(n=row_number(),
         x = (row_number() - 1) %% 12 + 1,
         y = ceiling(n  / 12)) %>%
  set_colnames(c("color","n","x","y")) %>% as_tibble()


#[1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF"
#[7] "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"

#[1] "#3B4992FF" "#EE0000FF" "#008B45FF" "#631879FF" "#008280FF" "#BB0021FF"
#[7] "#5F559BFF" "#A20056FF" "#808180FF" "#1B1919FF"

#[1] "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF"
#[7] "#FFDC91FF" "#EE4C97FF"

#[1] "#00468BFF" "#ED0000FF" "#42B540FF" "#0099B4FF" "#925E9FFF" "#FDAF91FF" "#AD002AFF" "#ADB6B6FF"
#[9] "#1B1919FF"

#[1] "#374E55FF" "#DF8F44FF" "#00A1D5FF" "#B24745FF" "#79AF97FF" "#6A6599FF" "#80796BFF"

#[1] "#0073C2FF" "#EFC000FF" "#868686FF" "#CD534CFF" "#7AA6DCFF" "#003C67FF" "#8F7700FF" "#3B3B3BFF"
#[9] "#A73030FF" "#4A6990FF"

#[1] "#FF6F00FF" "#C71000FF" "#008EA0FF" "#8A4198FF" "#5A9599FF" "#FF6348FF" "#84D7E1FF" "#FF95A8FF"
#[9] "#3D3B25FF" "#ADE2D0FF" "#1A5354FF" "#3F4041FF"

#[1] "#FF410DFF" "#6EE2FFFF" "#F7C530FF" "#95CC5EFF" "#D0DFE6FF" "#F79D1EFF" "#748AA6FF"

#[1] "#FED439FF" "#709AE1FF" "#8A9197FF" "#D2AF81FF" "#FD7446FF" "#D5E4A2FF" "#197EC0FF" "#F05C3BFF" "#46732EFF"
#[10] "#71D0F5FF" "#370335FF" "#075149FF" "#C80813FF" "#91331FFF" "#1A9993FF" "#FD8CC1FF"

#[1] "#CC0C00FF" "#5C88DAFF" "#84BD00FF" "#FFCD00FF" "#7C878EFF" "#00B5E2FF" "#00AF66FF"

#[1] "#2E2A2BFF" "#CF4E9CFF" "#8C57A2FF" "#358DB9FF" "#82581FFF" "#2F509EFF" "#E5614CFF" "#97A1A7FF" "#3DA873FF"
#[10] "#DC9445FF"

#[1] "#FAFD7CFF" "#82491EFF" "#24325FFF" "#B7E4F9FF" "#FB6467FF" "#526E2DFF" "#E762D7FF" "#E89242FF" "#FAE48BFF"
#[10] "#A6EEE6FF" "#917C5DFF" "#69C8ECFF"

#[1] "#C0392BFF" "#D35400FF" "#F39C12FF" "#27AE60FF" "#16A085FF" "#2980B9FF" "#8E44ADFF" "#2C3E50FF" "#7F8C8DFF"
#[10] "#BDC3C7FF"

#[1] "#D51317FF" "#F39200FF" "#EFD500FF" "#95C11FFF" "#007B3DFF" "#31B7BCFF" "#0094CDFF" "#164194FF" "#6F286AFF"
#[10] "#706F6FFF"




library(ggsci)

p <- ggplot(data=AI_GSE16561_data, 
            aes(x = Genes, y = Importance, fill = Genes)) +
  geom_bar(stat="identity") +
  #利用ggsci包中的调色板
  scale_fill_manual(values = df$color) +
  labs(x = "",title = 'GSE16561') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1)) +
  facet_wrap(~Method, scales = "free_y",ncol = 1)

ggpreview(p,width = 5,height = 6)
ggsave('AI_GSE16561_data.pdf',p,width = 5,height = 6)
ggsave('AI_GSE16561_data.tiff',p,width = 5,height = 6,dpi = 300)









##2.Immunoinfiltration analysis----

#创建新文件夹
dir.create("./Output/1_DEGs_WGCNA_GSEA/TME")

#更改路径至新文件夹
setwd("./Output/1_DEGs_WGCNA_GSEA/TME")


#**针对SD数据进行免疫浸润分析**


load('1_DEGs_WGCNA_GSEA/01_Exp_SD.Rdata')

IDs = c('GSE208668','GSE56931','GSE98582')

SD_Data_lists = list(Exp_GSE208668_Array_filter,
                     Exp_GSE56931_Array_filter,
                     Exp_GSE98582_Array_filter)

SD_Group_lists = list(Pd_GSE208668_filter %>% rownames_to_column('ID') %>% dplyr::select(ID,Group),
                      Pd_GSE56931 %>% dplyr::select(-ID) %>% rownames_to_column('ID') %>% dplyr::select(ID,Group),
                      Pd_GSE98582 %>% dplyr::select(-ID) %>% rownames_to_column('ID') %>% dplyr::select(ID,Group))


for (i in 1:length(SD_Data_lists)) {
  
  dat <- apply(SD_Data_lists[[i]], 2, as.numeric)  #因为此时的矩阵属于整数型，运行报错，因此需要先转化围殴数值型
  
  rownames(dat) <- rownames(SD_Data_lists[[i]])
  
  eset_stad = dat
  
  eset_stad = eset_stad[!duplicated(rownames(eset_stad)),]
  
  cibersort<-deconvo_tme(eset = eset_stad, method = "cibersort", arrays = T, perm = 5 )
  
  epic<-deconvo_tme(eset = eset_stad, method = "epic", arrays = T)
  
  mcp<-deconvo_tme(eset = eset_stad, method = "mcpcounter")
  
  xcell<-deconvo_tme(eset = eset_stad, method = "xcell",arrays = T)
  
  estimate<-deconvo_tme(eset = eset_stad, method = "estimate")
  
  timer<-deconvo_tme(eset = eset_stad, method = "timer", group_list = rep("stad",dim(eset_stad)[2]))
  
  quantiseq<-deconvo_tme(eset = eset_stad, tumor = TRUE, arrays = T, scale_mrna = TRUE, method = "quantiseq")
  
  ips<-deconvo_tme(eset = eset_stad, method = "ips", plot= FALSE)
  
  
  #**Combination of above deconvolution results**
  tme_combine<-cibersort %>% 
    inner_join(.,mcp,by       = "ID") %>% 
    inner_join(.,xcell,by     = "ID") %>%
    inner_join(.,epic,by      = "ID") %>% 
    inner_join(.,estimate,by  = "ID") %>% 
    inner_join(.,timer,by     = "ID") %>% 
    inner_join(.,quantiseq,by = "ID") %>% 
    inner_join(.,ips,by       = "ID")
  
  write.csv(tme_combine,file = paste0("TME/", 'SD_',IDs[[i]],"_TME_combine_result.csv"))
  
  
  #获取分组信息
  n <- SD_Group_lists[[i]]
  
  #定义不同方法的计算结果
  CIBERSORT <- (colnames(tme_combine)[grepl(colnames(tme_combine),pattern = 'CIBERSORT')])[1:22]
  EPIC <- colnames(tme_combine)[grepl(colnames(tme_combine),pattern = 'EPIC')]
  MCPcounter <- colnames(tme_combine)[grepl(colnames(tme_combine),pattern = 'MCPcounter')]
  xCell <- colnames(tme_combine)[grepl(colnames(tme_combine),pattern = 'xCell')] 
  Estimate <- colnames(tme_combine)[grepl(colnames(tme_combine),pattern = 'estimate')] 
  TIMER <- colnames(tme_combine)[grepl(colnames(tme_combine),pattern = 'TIMER')]
  quantiseq <- colnames(tme_combine)[grepl(colnames(tme_combine),pattern = 'quantiseq')]
  IPS <- (colnames(tme_combine)[grepl(colnames(tme_combine),pattern = 'IPS')])[1:6]
  
  Methods <- list(CIBERSORT,EPIC,MCPcounter,xCell,Estimate,TIMER,quantiseq,IPS)
  
  names(Methods) <- c('CIBERSORT','EPIC','MCPcounter','xCell','Estimate','TIMER','quantiseq','IPS')
  
  for (k in 1:length(Methods)) {
    
    res <- iobr_cor_plot(pdata_group           = n,
                         id1                   = "ID",
                         feature_data          = tme_combine,
                         id2                   = "ID",
                         target                = NULL,
                         group                 = "Group",#定义组别
                         is_target_continuous  = FALSE,
                         padj_cutoff           = 1,
                         index                 = 1,
                         category              = "signature",
                         signature_group       = Methods[k],
                         ProjectID             = names(Methods)[k],
                         palette_box           = "paired1",
                         palette_corplot       = "pheatmap",
                         palette_heatmap       = 2,
                         feature_limit         = 26,
                         character_limit       = 30,
                         show_heatmap_col_name = FALSE,
                         show_col              = FALSE,
                         show_plot             = TRUE,
                         path                  = paste0(IDs[i],'_',paste0(names(Methods)[k], "_SD_signatures")))
    
    write.csv(res,file = paste0("TME/",  paste0(IDs[i],'_',paste0(names(Methods)[k], "_SD_signatures")),"/result.csv"))
  }
}






##3.Characteristic gene set enrichment analysis----

rm(list = ls())

load('F:/01_Exp_SD.Rdata')

load('F:/01_Exp_Stroke.Rdata')


library(GSVA)
library(GSEABase)


#自定义函数

Signature_limma_dif <- function(ID = ID, Type = type, Matrix = matrix, Names_Matrix = Names_Matrix, Names_gemts = Names_gemts, 
                                Idtype = idtype, Org = org, gmts = gmts, Group = group, contrasts = contrasts){
  
  library(limma)
  
  if (Type == "RNA-seq") {
    
    eset_stad <- count2tpm(countMat = Matrix, idType = Idtype, org = Org, source = "local" )#目前支持hsa和mmu
    
    print("count2tpm finished")
    
    es_max <- lapply(gmts, function(gmtfile){ 
      
      geneset <- getGmt(file.path(here('Data/Subtype Features/'),gmtfile))  
      
      es_max <- gsva(as.matrix(eset_stad), 
                     
                     geneset,
                     
                     min.sz=5,max.sz=5000,
                     
                     mx.diff=FALSE, verbose=FALSE)
      
      return(es_max)})
    
    print("signature score finished") 
    
    es_dif <- lapply(es_max, function(es){ 
      
      es = es$es.obs
      
      design <- model.matrix(~0 + Group)
      
      colnames(design) = levels(factor(Group))
      
      rownames(design) = colnames(es)
      
      compare <- makeContrasts(contrasts = contrasts, levels=design)
      
      fit <- lmFit(es, design)
      
      fit2 <- contrasts.fit(fit, compare)
      
      fit3 <- eBayes(fit2)
      
      Diff <- topTable(fit3, coef=1, number=200)
      
      return(Diff)})
    
    save(es_max, es_dif, file =  here(paste0('Data/', ID ), paste0(Names_Matrix, '_', 'es_max_dif.Rdata')))
    
    library(openxlsx)
    
    wb_dif <- createWorkbook()
    
    for (i in 1:length(es_dif)) {
      
      # 创建一个新的sheet，并以元素的名称命名
      addWorksheet(wb_dif, sheetName = Names_gemts[i])
      
      # 将元素写入当前sheet
      writeData(wb_dif, sheet = i, x = es_dif[[i]], rowNames = T)
    }
    
    # 保存Excel文件
    saveWorkbook(wb_dif, file = here(paste0('Output/', ID ), paste0(Names_Matrix, '_', "es_dif_output.xlsx")), overwrite = TRUE)
    
    print("signature Dif finished") 
    
  }
  
  if (Type == "Array") {
    
    eset_stad <- Matrix#目前支持hsa和mmu
    
    print("Array finished")
    
    es_max <- lapply(gmts, function(gmtfile){ 
      
      geneset <- getGmt(file.path(here('Data/Subtype Features/'),gmtfile))  
      
      es_max <- gsva(as.matrix(eset_stad), 
                     
                     geneset,
                     
                     min.sz=5,max.sz=5000,
                     
                     mx.diff=FALSE, verbose=FALSE)
      
      return(es_max)})
    
    print("signature score finished") 
    
    es_dif <- lapply(es_max, function(es){ 
      
      es = es$es.obs
      
      design <- model.matrix(~0 + Group)
      
      colnames(design) = levels(factor(Group))
      
      rownames(design) = colnames(es)
      
      compare <- makeContrasts(contrasts = contrasts, levels=design)
      
      fit <- lmFit(es, design)
      
      fit2 <- contrasts.fit(fit, compare)
      
      fit3 <- eBayes(fit2)
      
      Diff <- topTable(fit3, coef=1, number=200)
      
      return(Diff)})
    
    save(es_max, es_dif, file =  here(paste0('Data/', ID ,'/'), paste0(Names_Matrix, '_', 'es_max_dif.Rdata')))
    
    library(openxlsx)
    
    wb_dif <- createWorkbook()
    
    for (i in 1:length(es_dif)) {
      
      # 创建一个新的sheet，并以元素的名称命名
      addWorksheet(wb_dif, sheetName = Names_gemts[i])
      
      # 将元素写入当前sheet
      writeData(wb_dif, sheet = i, x = es_dif[[i]], rowNames = T)
    }
    
    # 保存Excel文件
    saveWorkbook(wb_dif, file = here(paste0('Output/', ID ), paste0(Names_Matrix, '_', "es_dif_output.xlsx")), overwrite = TRUE)
    
    print("signature Dif finished") 
    
  }
}



gmts=list.files('',pattern = 'gmt')

gmts
# [1] "c2.cp.kegg.v2022.1.Hs.symbols.gmt"         "c2.cp.reactome.v2022.1.Hs.symbols.gmt"    
# [3] "c2.cp.wikipathways.v2022.1.Hs.symbols.gmt" "c5.go.bp.v2022.1.Hs.symbols.gmt"          
# [5] "Cancer_Cell_2021_Features.symbols.gmt"     "h.all.v2022.1.Hs.symbols.gmt"             
# [7] "IOBR_Features.symbols.gmt"   

names_gemts <- c("kegg","reactome",
                 "wikipathways","go.bp",
                 "Cancer_Cell","h.all","IOBR_Features")


#补充分组
GSE208668_Group <- Pd_GSE208668_filter$Group

GSE16561_Group <- Pd_GSE16561$Group

Signature_limma_dif(ID = '1_DEGs_WGCNA_GSEA', Type = 'Array', Matrix = Exp_GSE208668_Array_filter,
                    Names_gemts = names_gemts,Names_Matrix = 'Exp_GSE208668',
                    Idtype = 'Symbol', Org = 'hsa', gmts = gmts, Group = GSE208668_Group, contrasts = 'Control-SD')

Signature_limma_dif(ID = '1_DEGs_WGCNA_GSEA', Type = 'Array', Matrix = Exp_GSE16561_Array_filter,
                    Names_gemts = names_gemts,Names_Matrix = 'Exp_GSE16561',
                    Idtype = 'Symbol', Org = 'hsa', gmts = gmts, Group = GSE16561_Group, contrasts = 'Control-Stroke')




#挑选特征通路进行展示

load(here('/Exp_GSE208668_es_max_dif.Rdata'))

GSE208668_data1 = es_max[[5]]$es.obs %>% 
  as.data.frame()

GSE208668_data2 = es_max[[6]]$es.obs %>% 
  as.data.frame()

annotation_col = data.frame(Group = ifelse(Pd_GSE208668_filter$Group=='SD','Control','SD'))
rownames(annotation_col) = rownames(Pd_GSE208668_filter)

ann_colors = list(Group = c(Control = "#20854EFF", SD = "#E18727FF"))

pdf('Heatmap_GSE208668_Enrich1.pdf',width = 10,height = 5)
pheatmap(GSE208668_data1, 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cellwidth = 12, cellheight = 12,
         show_colnames = F,border_color = 'white',
         cluster_rows = F,
         cluster_cols = F)
dev.off()


pdf('Heatmap_GSE208668_Enrich2.pdf',width = 14,height = 10)
pheatmap(GSE208668_data2, 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cellwidth = 12, cellheight = 12,
         show_colnames = F,border_color = 'white',
         cluster_rows = F,
         cluster_cols = F)
dev.off()





load(here('Data/1_DEGs_WGCNA_GSEA/Exp_GSE16561_es_max_dif.Rdata'))

GSE16561_data1 = es_max[[5]]$es.obs %>% 
  as.matrix()

GSE16561_data2 = es_max[[6]]$es.obs %>% 
  as.matrix()

annotation_col = data.frame(Group = ifelse(Pd_GSE16561$Group=='Control','Control','Stroke') %>%
                              factor(levels = c('Control','Stroke')))
rownames(annotation_col) = rownames(Pd_GSE16561)

ann_colors = list(Group = c(Control = "#008EA0FF", Stroke = "#C71000FF"))

pdf('Heatmap_GSE16561_Enrich1.pdf',width = 10,height = 5)
pheatmap(GSE16561_data1, 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cellwidth = 12, cellheight = 12,
         show_colnames = F,border_color = 'white',
         cluster_rows = F,
         cluster_cols = F)
dev.off()


pdf('Heatmap_GSE16561_Enrich2.pdf',width = 16,height = 12)
pheatmap(GSE16561_data2, 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cellwidth = 12, cellheight = 12,
         show_colnames = F,border_color = 'white',
         cluster_rows = F,
         cluster_cols = F)
dev.off()





##4.ROC and expression analysis----

Genes
# [1] "PUF60"     "IMP4"      "NDUFB8"    "REXO4"     "DGUOK"     "CECR5"     "ARL2"     
# [8] "ZC3HC1"    "KLHL22"    "ADSL"      "HDAC1"     "ETV3"      "ORAOV1"    "PUS1"     
# [15] "SRPRB"     "P2RY11"    "CCDC12"    "RSL1D1"    "LTA"       "NSMCE1"    "SLC41A3"  
# [22] "CD81"      "PPP1R13B"  "CCDC92"    "WDR54"     "FNTA"      "RCC2"      "RAB11FIP3"
# [29] "ACACB"     "ZNF671"    "GRAP"      "FN3KRP"    "SCAMP3"    "MLLT6"     "AP3M2"    
# [36] "CCDC84"    "BAG3" 


library(tidyplots)


GSE208668_datExpr <- Exp_GSE208668_Array_filter %>% 
  as.data.frame() %>% 
  filter(row.names(.) %in% Genes) %>%
  rownames_to_column('Genes') %>% 
  pivot_longer(cols = -1,names_to = 'Sample',values_to = 'Expression') %>%
  mutate(Group = ifelse(Sample %in% rownames(Pd_GSE208668_filter[Pd_GSE208668_filter$Group=='SD',]),'SD','Control') %>% 
           factor(levels = c('Control','SD')))

p = GSE208668_datExpr %>% 
  tidyplot(x = Group, y = Expression, color = Group) %>% 
  add_mean_dash() %>% 
  add_sem_errorbar() %>% 
  add_data_points() %>%  
  remove_legend() %>% 
  adjust_x_axis_title('') %>%
  adjust_y_axis_title('Scale_Expression') %>% 
  adjust_colors(new_colors = c("#20854EFF","#E18727FF")) %>%
  add_test_asterisks(ref.group = 'Control',hide_info = TRUE) %>% 
  split_plot(by = Genes,nrow = 5) 
ggpreview(p, width = 14, height = 8)
#ggsave(here(','GSE208668_Genes_boxplot.pdf'),p,width = 14, height = 8)
#ggsave(here(','GSE208668_Genes_boxplot.tiff'),p,width = 14, height = 8,dpi = 300)




GSE16561_datExpr <- Exp_GSE16561_Array_filter %>% 
  as.data.frame() %>% 
  filter(row.names(.) %in% Genes) %>%
  rownames_to_column('Genes') %>% 
  pivot_longer(cols = -1,names_to = 'Sample',values_to = 'Expression') %>%
  mutate(Group = ifelse(Sample %in% rownames(Pd_GSE16561[Pd_GSE16561$Group=='Stroke',]),'Stroke','Control') %>% 
           factor(levels = c('Control','Stroke')))

p = GSE16561_datExpr %>% 
  tidyplot(x = Group, y = Expression, color = Group) %>% 
  add_mean_dash() %>% 
  add_sem_errorbar() %>% 
  add_data_points() %>%  
  remove_legend() %>% 
  adjust_x_axis_title('') %>%
  adjust_colors(new_colors = c("#008EA0FF","#C71000FF")) %>%
  adjust_y_axis_title('Scale_Expression') %>% 
  add_test_asterisks(ref.group = 'Control',hide_info = TRUE) %>% 
  split_plot(by = Genes,nrow = 5) 
ggpreview(p, width = 14, height = 8)
#ggsave(here(','GSE16561_Genes_boxplot.pdf'),p,width = 14, height = 8)
#ggsave(here(','GSE16561_Genes_boxplot.tiff'),p,width = 14, height = 8,dpi = 300)





library(ggplotify)
library(cowplot)
library(pROC)

data = Exp_GSE208668_Array_filter %>% 
  as.data.frame() %>% 
  filter(row.names(.) %in% Genes) %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = 'ID') %>%
  mutate(Group = Pd_GSE208668_filter$Group)

#构建循环，针对每一个基因进行ROC，然后输出AUC值组成新的数据框

roc_data = data.frame()

for (i in 1:nrow(data)) {
  
  fit = roc(data$Group, data[,i+1])
  
  roc_data = rbind(roc_data,data.frame(Genes = colnames(data)[i+1],AUC = fit$auc))
}

fit = roc(data$Group, data$gene_name)

plot.roc(fit,
         axes=T, ## 是否显示xy轴
         legacy.axes=T,
         main='ARL2', ## Title
         col= "#91D1C2FF", ## 曲线颜色
         lty=1, ## 曲线形状
         lwd=3, ## 曲线粗细
         identity=T, ## 是否显示对角线
         identity.col="grey60", ## 对角线颜色
         identity.lty=2, ## 对角线形状
         identity.lwd=2, ## 对角线粗细
         print.thres=F, ## 是否输出cut-off值
         print.thres.pch=20, ## cut-off点的形状
         print.thres.col="red", ## cut-off点和文本的颜色
         print.thres.cex=1.2, 
         print.auc=T, ## 是否显示AUC
         print.auc.pattern=paste0("AUC = ",round(fit$auc,2)), ## 展示AUC的格式
         auc.polygon.border="#91D1C2FF",
         print.auc.x=0.5, ## AUC值的X位置
         print.auc.y=0.5, ## AUC值的Y位置
         print.auc.cex=1, ## AUC值的放大倍数
         print.auc.col='black', ## ACU值的颜色
         auc.polygon=TRUE, ## 是否将ROC下面积上色
         auc.polygon.col="#91D1C2FF", 
         max.auc.polygon=F,
         max.auc.polygon.col='WhiteSmoke',
         max.auc.polygon.lty=1
)

p <-recordPlot()#函数截取画布内容

p<- as.ggplot(ggdraw(p))#draw绘制截取内容并将其转化为ggplot2类型

ggsave(file = "Gene_ROC.pdf",p,width = 4,height = 4)




data = Exp_GSE16561_Array_filter %>% 
  as.data.frame() %>% 
  filter(row.names(.) %in% Genes) %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = 'ID') %>%
  mutate(Group = Pd_GSE16561$Group)


roc_data = data.frame()

for (i in 1:nrow(data)) {
  
  fit = roc(data$Group, data[,i+1])
  
  roc_data = rbind(roc_data,data.frame(Genes = colnames(data)[i+1],AUC = fit$auc))
}

fit = roc(data$Group, data$gene_name)

plot.roc(fit,
         axes=T, ## 是否显示xy轴
         legacy.axes=T,
         main='ARL2', ## Title
         col= "#FFDC91FF", ## 曲线颜色
         lty=1, ## 曲线形状
         lwd=3, ## 曲线粗细
         identity=T, ## 是否显示对角线
         identity.col="grey60", ## 对角线颜色
         identity.lty=2, ## 对角线形状
         identity.lwd=2, ## 对角线粗细
         print.thres=F, ## 是否输出cut-off值
         print.thres.pch=20, ## cut-off点的形状
         print.thres.col="red", ## cut-off点和文本的颜色
         print.thres.cex=1.2, 
         print.auc=T, ## 是否显示AUC
         print.auc.pattern=paste0("AUC = ",round(fit$auc,2)), ## 展示AUC的格式
         auc.polygon.border="#FFDC91FF",
         print.auc.x=0.5, ## AUC值的X位置
         print.auc.y=0.5, ## AUC值的Y位置
         print.auc.cex=1, ## AUC值的放大倍数
         print.auc.col='black', ## ACU值的颜色
         auc.polygon=TRUE, ## 是否将ROC下面积上色
         auc.polygon.col="#FFDC91FF", 
         max.auc.polygon=F,
         max.auc.polygon.col='WhiteSmoke',
         max.auc.polygon.lty=1
)

p <-recordPlot()#函数截取画布内容

p<- as.ggplot(ggdraw(p))#draw绘制截取内容并将其转化为ggplot2类型

ggsave(file = "Gene_ROC.pdf",p,width = 4,height = 4)










data = Exp_GSE208668_Array_filter %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% c('MICA','GPR68','MLKL','PTGDR',"OMD","FKBP5",'CPM',
                            "ACACB","ZNF671","GRAP","FN3KRP","SCAMP3","MLLT6","AP3M2")) %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = 'ID') %>%
  mutate(Group = Pd_GSE208668_filter$Group) %>% 
  pivot_longer(cols = -c(ID,Group),names_to = 'Genes',values_to = 'Expression')


Dat = data %>% 
  pivot_wider(names_from = Genes,values_from = Expression)


#绘制ROC曲线

fit = roc(Dat$Group, Dat$gene_name)

plot.roc(fit,
         axes=T, ## 是否显示xy轴
         legacy.axes=T,
         main='ARL2', ## Title
         col= "#95CC5EFF", ## 曲线颜色
         lty=1, ## 曲线形状
         lwd=3, ## 曲线粗细
         identity=T, ## 是否显示对角线
         identity.col="grey60", ## 对角线颜色
         identity.lty=2, ## 对角线形状
         identity.lwd=2, ## 对角线粗细
         print.thres=F, ## 是否输出cut-off值
         print.thres.pch=20, ## cut-off点的形状
         print.thres.col="red", ## cut-off点和文本的颜色
         print.thres.cex=1.2, 
         print.auc=T, ## 是否显示AUC
         print.auc.pattern=paste0("AUC = ",round(fit$auc,2)), ## 展示AUC的格式
         auc.polygon.border="#95CC5EFF",
         print.auc.x=0.5, ## AUC值的X位置
         print.auc.y=0.5, ## AUC值的Y位置
         print.auc.cex=1, ## AUC值的放大倍数
         print.auc.col='black', ## ACU值的颜色
         auc.polygon=TRUE, ## 是否将ROC下面积上色
         auc.polygon.col="#95CC5EFF", 
         max.auc.polygon=F,
         max.auc.polygon.col='WhiteSmoke',
         max.auc.polygon.lty=1
)

p <-recordPlot()#函数截取画布内容

p <- as.ggplot(ggdraw(p))#draw绘制截取内容并将其转化为ggplot2类型

ggsave(file = "Gene_ROC.pdf",p,width = 4,height = 4)


#绘制基因表达差异箱线图
p = data %>% 
  filter(Genes == 'CPM') %>%
  tidyplot(x = Group, y = Expression, color = Group) %>% 
  add_mean_dash() %>% 
  add_sem_errorbar() %>% 
  add_data_points() %>%  
  adjust_x_axis_title('') %>%
  adjust_colors(new_colors = c("#20854EFF","#E18727FF")) %>%
  adjust_y_axis_title('Scale_Expression') %>% adjust_title('GSE22255') %>%
  add_test_asterisks(ref.group = 'Control',hide_info = TRUE)

ggpreview(p, width = 4, height = 4)
ggsave()






fit = roc()

plot.roc(fit,
         axes=T, ## 是否显示xy轴
         legacy.axes=T,
         main='ARL2', ## Title
         col= "#F39B7FFF", ## 曲线颜色
         lty=1, ## 曲线形状
         lwd=3, ## 曲线粗细
         identity=T, ## 是否显示对角线
         identity.col="grey60", ## 对角线颜色
         identity.lty=2, ## 对角线形状
         identity.lwd=2, ## 对角线粗细
         print.thres=F, ## 是否输出cut-off值
         print.thres.pch=20, ## cut-off点的形状
         print.thres.col="red", ## cut-off点和文本的颜色
         print.thres.cex=1.2, 
         print.auc=T, ## 是否显示AUC
         print.auc.pattern=paste0("AUC = ",round(fit$auc,2)), ## 展示AUC的格式
         auc.polygon.border="#F39B7FFF",
         print.auc.x=0.5, ## AUC值的X位置
         print.auc.y=0.5, ## AUC值的Y位置
         print.auc.cex=1, ## AUC值的放大倍数
         print.auc.col='black', ## ACU值的颜色
         auc.polygon=TRUE, ## 是否将ROC下面积上色
         auc.polygon.col="#F39B7FFF", 
         max.auc.polygon=F,
         max.auc.polygon.col='WhiteSmoke',
         max.auc.polygon.lty=1
)

p <-recordPlot()#函数截取画布内容

p <- as.ggplot(ggdraw(p))#draw绘制截取内容并将其转化为ggplot2类型

ggsave(file = "GSE98566_Gene_ROC.pdf",p,width = 4,height = 4)


#绘制基因表达差异箱线图
p = data %>% 
  filter(Genes == '') %>%
  tidyplot(x = Group, y = Expression, color = Group) %>% 
  add_mean_dash() %>% 
  add_sem_errorbar() %>% 
  add_data_points() %>%  
  adjust_x_axis_title('') %>%
  adjust_colors(new_colors = c("#008EA0FF","#C71000FF")) %>%
  adjust_y_axis_title('Scale_Expression') %>% adjust_title('GSE98566') %>%
  add_test_asterisks(ref.group = 'Control',hide_info = TRUE)

ggpreview(p, width = 4, height = 4)
#ggsave(here(','GSE98566_Genes_boxplot.pdf'),p,width = 4, height = 4)
#ggsave(here(','GSE98566_Genes_boxplot.tiff'),p,width = 4, height = 4,dpi = 300)






