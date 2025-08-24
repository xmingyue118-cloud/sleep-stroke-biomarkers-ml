



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


##############################################################
################SD data sets######################################
##############################################################


#GSE208668 data pre-processing----

library(vroom)
library(tidyverse)

##1.Extract the expression matrix----

library(GEOquery)
Eset <- getGEO("GSE208668",destdir = '//GSE208668_Rawdata')
expr_counts <- Eset$GSE208668_series_matrix.txt.gz@assayData$exprs
Pd <- Eset$GSE208668_series_matrix.txt.gz@phenoData@data

colnames(Pd)
# [1] "title"                     "geo_accession"             "status"                   
# [4] "submission_date"           "last_update_date"          "type"                     
# [7] "channel_count"             "source_name_ch1"           "organism_ch1"             
# [10] "characteristics_ch1"       "characteristics_ch1.1"     "characteristics_ch1.2"    
# [13] "characteristics_ch1.3"     "characteristics_ch1.4"     "characteristics_ch1.5"    
# [16] "characteristics_ch1.6"     "characteristics_ch1.7"     "characteristics_ch1.8"    
# [19] "characteristics_ch1.9"     "molecule_ch1"              "extract_protocol_ch1"     
# [22] "label_ch1"                 "label_protocol_ch1"        "taxid_ch1"                
# [25] "hyb_protocol"              "scan_protocol"             "description"              
# [28] "data_processing"           "platform_id"               "contact_name"             
# [31] "contact_email"             "contact_department"        "contact_institute"        
# [34] "contact_address"           "contact_city"              "contact_state"            
# [37] "contact_zip/postal_code"   "contact_country"           "supplementary_file"       
# [40] "data_row_count"            "age:ch1"                   "bdi:ch1"                  
# [43] "bdins:ch1"                 "bmi:ch1"                   "comorbidity:ch1"          
# [46] "education (years):ch1"     "gender:ch1"                "history of depression:ch1"
# [49] "SD:ch1"              "race:ch1"     

#View sample distribution
table(Pd$`insomnia:ch1`)
# no yes 
# 25  17 

library(limma)

expr_counts <- normalizeBetweenArrays(expr_counts)

boxplot(expr_counts[,1:10])


#No probe conversion is required





##2.Organize clinical information----


colnames(Pd)
# [1] "title"                     "geo_accession"             "status"                   
# [4] "submission_date"           "last_update_date"          "type"                     
# [7] "channel_count"             "source_name_ch1"           "organism_ch1"             
# [10] "characteristics_ch1"       "characteristics_ch1.1"     "characteristics_ch1.2"    
# [13] "characteristics_ch1.3"     "characteristics_ch1.4"     "characteristics_ch1.5"    
# [16] "characteristics_ch1.6"     "characteristics_ch1.7"     "characteristics_ch1.8"    
# [19] "characteristics_ch1.9"     "molecule_ch1"              "extract_protocol_ch1"     
# [22] "label_ch1"                 "label_protocol_ch1"        "taxid_ch1"                
# [25] "hyb_protocol"              "scan_protocol"             "description"              
# [28] "data_processing"           "platform_id"               "contact_name"             
# [31] "contact_email"             "contact_department"        "contact_institute"        
# [34] "contact_address"           "contact_city"              "contact_state"            
# [37] "contact_zip/postal_code"   "contact_country"           "supplementary_file"       
# [40] "data_row_count"            "age:ch1"                   "bdi:ch1"                  
# [43] "bdins:ch1"                 "bmi:ch1"                   "comorbidity:ch1"          
# [46] "education (years):ch1"     "gender:ch1"                "history of depression:ch1"
# [49] "SD:ch1"              "race:ch1"     


Pd <- Pd %>% 
  as.data.frame() %>%
  mutate(Group = ifelse(`insomnia:ch1` == 'no', 'Control', 'SD') %>% factor(levels = c('Control','SD'))) %>% 
  dplyr::select(c("age:ch1","bdi:ch1","bmi:ch1","comorbidity:ch1",          
                  "education (years):ch1","gender:ch1","history of depression:ch1",
                  "race:ch1","Group")) %>% 
  dplyr::rename("Age" = "age:ch1",
                "BDI" = "bdi:ch1",
                "BMI" = "bmi:ch1",
                "Comorbidity" = "comorbidity:ch1",
                "Education" = "education (years):ch1",
                "Gender" = "gender:ch1",
                "Depression" = "history of depression:ch1",
                "Race" = "race:ch1")

Pd[1:5,]
# Age BDI         BMI Comorbidity Education Gender Depression      Race Group
# GSM6360934  65  13 21.49923325 0.638977647        16 female        yes     white    SD
# GSM6360935  75   7 26.41070366  0.95846647        16   male         no     white    SD
# GSM6360936  77   4 31.28330994  1.91693294        15 female        yes     white    SD
# GSM6360937  64   7  25.7443676           0        16 female         no non-white    SD
# GSM6360938  60   0 31.59882355           0        16   male        yes     white    SD

Exp_GSE208668_Array <- expr_counts

Pd_GSE208668 <- Pd



##3.PCA----


###**PCA**

#Set up grouping
Group <- Pd_GSE208668$Group %>% factor(levels = c('Control','SD'))

#Start by filtering the gene to ensure that it is expressed in at least half of the sample

table(rowSums(Exp_GSE208668_Array>=1)>=round(ncol(Exp_GSE208668_Array)/2))
# TRUE 
# 33209  

Exp_GSE208668_Array_filter <- Exp_GSE208668_Array[rowSums(Exp_GSE208668_Array>=1)>=round(ncol(Exp_GSE208668_Array)/2),]
dim(Exp_GSE208668_Array_filter)
# [1] 33210    42



#Principal component analysis
set.seed(10086)
dat = as.matrix(t(Exp_GSE208668_Array_filter))

dat = dat[,colSums(is.na(dat))<nrow(dat)]
pca<-prcomp(dat,scale. = T)
pca.results <- data.frame(pca$x)[,1:2]
pca.results$Group<- Group
pca.results[1:4,1:3]
# PC1      PC2    Group
# GSM6360934 -62.97908 57.86842 SD
# GSM6360935 -61.84146 56.78580 SD
# GSM6360936 -62.56079 53.31404 SD
# GSM6360937 -60.91480 63.11947 SD


centroid <- aggregate(cbind(PC1,PC2) ~ Group,
                      data = pca.results,
                      FUN = mean)

levels(centroid$Group)
# [1] "Control"  "SD"

centroid$PC1
# [1]  43.05739 -63.31970

centroid$PC2
# [1] -38.90453  57.21254



pca.results1 <- dplyr::left_join(pca.results, centroid, by = "Group", 
                                 suffix = c("",".cen"))

cols <- c('#37629C','#F18625')



pca_eig1 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[1,2],2)
pca_eig2 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[2,2],2)



PCA <- ggplot(pca.results1,aes(x=PC1,y=PC2))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(aes(color=Group),size = 3)+
  geom_segment(aes(xend=PC1.cen,yend=PC2.cen,color=Group),
               show.legend = F)+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', size = 1, fill = 'transparent'), 
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.5,"cm"),
        legend.box.spacing = unit(0,"cm"))+
  labs(x =  paste('PCA_1:', pca_eig1, '%'), y = paste('PCA_2:', pca_eig2, '%'), color = '', title = 'GSE208668')+#将 PCA 轴贡献度添加到坐标轴标题中
  scale_colour_manual(values = cols)+
  geom_label(data = centroid, 
             aes(label = Group, fill = Group), size = 4, 
             show.legend = FALSE,
             color="white")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "top")


ggpreview(PCA,width = 3,height = 3.5)
ggsave(PCA,filename = here('','PCA_GSE208668_All.pdf'),width = 3,height = 3.5)
ggsave(PCA,filename = here('','PCA_GSE208668_All.tiff'),width = 3,height = 3.5,dpi = 300)





#According to the PCA data, some samples were deleted
Exp_GSE208668_Array_filter <- Exp_GSE208668_Array_filter[,!colnames(Exp_GSE208668_Array_filter)%in%rownames(pca.results[pca.results$PC1>0 & pca.results$Group=='Control',])] %>% 
  as.data.frame() %>% 
  na.omit()

dim(Exp_GSE208668_Array_filter)
# [1] 33209    33

Pd_GSE208668_filter <- Pd_GSE208668[match(colnames(Exp_GSE208668_Array_filter),rownames(Pd_GSE208668)),]




set.seed(10086)
dat = as.matrix(t(Exp_GSE208668_Array_filter))

dat = dat[,colSums(is.na(dat))<nrow(dat)]
pca<-prcomp(dat,scale. = T)
pca.results <- data.frame(pca$x)[,1:2]
pca.results$Group<- Pd_GSE208668_filter$Group
pca.results[1:4,1:3]
# PC1       PC2 Group
# GSM6360934 -72.36380   3.66478    SD
# GSM6360935 -70.86557 104.14033    SD
# GSM6360936 -65.94930 110.43308    SD
# GSM6360937 -85.86780  77.00190    SD


centroid <- aggregate(cbind(PC1,PC2) ~ Group,
                      data = pca.results,
                      FUN = mean)

levels(centroid$Group)
# [1] "Control"  "SD"

centroid$PC1
# [1]  80.76434 -76.01350

centroid$PC2




pca.results1 <- dplyr::left_join(pca.results, centroid, by = "Group", 
                                 suffix = c("",".cen"))

cols <- c('#37629C','#F18625')



pca_eig1 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[1,2],2)
pca_eig2 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[2,2],2)



PCA <- ggplot(pca.results1,aes(x=PC1,y=PC2))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(aes(color=Group),size = 3)+
  geom_segment(aes(xend=PC1.cen,yend=PC2.cen,color=Group),
               show.legend = F)+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', size = 1, fill = 'transparent'), 
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.5,"cm"),
        legend.box.spacing = unit(0,"cm"))+
  labs(x =  paste('PCA_1:', pca_eig1, '%'), y = paste('PCA_2:', pca_eig2, '%'), color = '', title = 'GSE208668')+#将 PCA 轴贡献度添加到坐标轴标题中
  scale_colour_manual(values = cols)+
  geom_label(data = centroid, 
             aes(label = Group, fill = Group), size = 4, 
             show.legend = FALSE,
             color="white")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "top")


ggpreview(PCA,width = 3,height = 3.5)
ggsave(PCA,filename = here('Output/1_DEGs_WGCNA_GSEA/','PCA_GSE208668.pdf'),width = 3,height = 3.5)
ggsave(PCA,filename = here('Output/1_DEGs_WGCNA_GSEA/','PCA_GSE208668.tiff'),width = 3,height = 3.5,dpi = 300)


#GSE98566 data pre-processing----

library(GEOquery)
Eset <- getGEO("GSE98566",destdir = 'C:/Users/XL/Documents/SD/Zhang_IS_SD/Data/Rawdata/GSE98566')
expr_counts <- Eset$GSE98566_series_matrix.txt.gz@assayData$exprs
Pd <- Eset$GSE98566_series_matrix.txt.gz@phenoData@data

gse_1 <- getGEO(filename = "GPL6244.annot.gz")
gse <- getGEO(filename = "GPL6244.soft.gz")

#View sample distribution
table(Pd$`subject group:ch1`)
#Control Sleep Deprived 
#71             92 


library(limma)

expr_counts <- normalizeBetweenArrays(expr_counts)

boxplot(expr_counts[,1:10])


#No probe conversion is required





##2.Organize clinical information----


colnames(Pd)
# [1] "title"                     "geo_accession"             "status"                   
# [4] "submission_date"           "last_update_date"          "type"                     
# [7] "channel_count"             "source_name_ch1"           "organism_ch1"             
# [10] "characteristics_ch1"       "characteristics_ch1.1"     "characteristics_ch1.2"    
# [13] "characteristics_ch1.3"     "characteristics_ch1.4"     "characteristics_ch1.5"    
# [16] "characteristics_ch1.6"     "characteristics_ch1.7"     "characteristics_ch1.8"    
# [19] "characteristics_ch1.9"     "molecule_ch1"              "extract_protocol_ch1"     
# [22] "label_ch1"                 "label_protocol_ch1"        "taxid_ch1"                
# [25] "hyb_protocol"              "scan_protocol"             "description"              
# [28] "data_processing"           "platform_id"               "contact_name"             
# [31] "contact_email"             "contact_department"        "contact_institute"        
# [34] "contact_address"           "contact_city"              "contact_state"            
# [37] "contact_zip/postal_code"   "contact_country"           "supplementary_file"       
# [40] "data_row_count"            "age:ch1"                   "bdi:ch1"                  
# [43] "bdins:ch1"                 "bmi:ch1"                   "comorbidity:ch1"          
# [46] "education (years):ch1"     "gender:ch1"                "history of depression:ch1"
# [49] "SD:ch1"              "race:ch1"     



Pd <- Pd %>% 
  as.data.frame() %>%mutate(Group = ifelse(`subject group:ch1` == 'Control', 'Control', 'SD') %>% factor(levels = c('Control','SD')))
# Age BDI         BMI Comorbidity Education Gender Depression      Race Group
# GSM6360934  65  13 21.49923325 0.638977647        16 female        yes     white    SD
# GSM6360935  75   7 26.41070366  0.95846647        16   male         no     white    SD
# GSM6360936  77   4 31.28330994  1.91693294        15 female        yes     white    SD
# GSM6360937  64   7  25.7443676           0        16 female         no non-white    SD
# GSM6360938  60   0 31.59882355           0        16   male        yes     white    SD

Exp_GSE98566_Array <- expr_counts

Pd_GSE98566 <- Pd



##3.PCA and correlation heatmaps----


###**PCA**

#Set up grouping
Group <- Pd_GSE98566$Group %>% factor(levels = c('Control','SD'))

#Start by filtering the gene to ensure that it is expressed in at least half of the sample

table(rowSums(Exp_GSE98566_Array>=1)>=round(ncol(Exp_GSE98566_Array)/2))
#TRUE 
#8515  

Exp_GSE98566_Array_filter <- Exp_GSE98566_Array[rowSums(Exp_GSE98566_Array>=1)>=round(ncol(Exp_GSE98566_Array)/2),]
dim(Exp_GSE98566_Array_filter)
#[1] 8515  163


#Principal component analysis
set.seed(10086)
dat = as.matrix(t(Exp_GSE98566_Array_filter))

dat = dat[,colSums(is.na(dat))<nrow(dat)]
pca<-prcomp(dat,scale. = T)
pca.results <- data.frame(pca$x)[,1:2]
pca.results$Group<- Group
pca.results[1:4,1:3]
#PC1       PC2   Group
#GSM2600155 -27.709953 -26.91261 Control
#GSM2600156 -19.642971  16.98514 Control
#GSM2600157   1.118866 -32.85220 Control
#GSM2600158 -21.039790 -15.04844 Control


centroid <- aggregate(cbind(PC1,PC2) ~ Group,
                      data = pca.results,
                      FUN = mean)

levels(centroid$Group)
# [1]  "Control" "SD" 

centroid$PC1
# [1]  -11.363055   8.769314

centroid$PC2



pca.results1 <- dplyr::left_join(pca.results, centroid, by = "Group", 
                                 suffix = c("",".cen"))

cols <- c('#37629C','#F18625')

library(FactoMineR)

pca_eig1 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[1,2],2)
pca_eig2 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[2,2],2)



PCA <- ggplot(pca.results1,aes(x=PC1,y=PC2))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(aes(color=Group),size = 3)+
  geom_segment(aes(xend=PC1.cen,yend=PC2.cen,color=Group),
               show.legend = F)+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', size = 1, fill = 'transparent'), 
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.5,"cm"),
        legend.box.spacing = unit(0,"cm"))+
  labs(x =  paste('PCA_1:', pca_eig1, '%'), y = paste('PCA_2:', pca_eig2, '%'), color = '', title = 'GSE98566')+#将 PCA 轴贡献度添加到坐标轴标题中
  scale_colour_manual(values = cols)+
  geom_label(data = centroid, 
             aes(label = Group, fill = Group), size = 4, 
             show.legend = FALSE,
             color="white")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "top")
PCA

ggpreview(PCA,width = 3,height = 3.5)
#ggsave(PCA,filename = 'PCA_GSE98566_All.pdf',width = 3,height = 3.5，dpi = 300)
#ggsave(PCA,filename = here('','PCA_GSE98566_All.tiff'),width = 3,height = 3.5,dpi = 300)





#According to the PCA data, some samples were deleted
Exp_GSE98566_Array_filter <- Exp_GSE98566_Array_filter[,!colnames(Exp_GSE98566_Array_filter)%in%rownames(pca.results[pca.results$PC1>0 & pca.results$Group=='Control',])] %>% 
  as.data.frame() %>% 
  na.omit()

dim(Exp_GSE98566_Array_filter)
# [1] 33209    33

Pd_GSE98566_filter <- Pd_GSE98566[match(colnames(Exp_GSE98566_Array_filter),rownames(Pd_GSE98566)),]




set.seed(10086)
dat = as.matrix(t(Exp_GSE98566_Array_filter))

dat = dat[,colSums(is.na(dat))<nrow(dat)]
pca<-prcomp(dat,scale. = T)
pca.results <- data.frame(pca$x)[,1:2]
pca.results$Group<- Pd_GSE98566_filter$Group
pca.results[1:4,1:3]
# PC1       PC2 Group
# GSM6360934 -72.36380   3.66478    SD
# GSM6360935 -70.86557 104.14033    SD
# GSM6360936 -65.94930 110.43308    SD
# GSM6360937 -85.86780  77.00190    SD


centroid <- aggregate(cbind(PC1,PC2) ~ Group,
                      data = pca.results,
                      FUN = mean)

levels(centroid$Group)
# [1] "Control"  "SD"

centroid$PC1
# [1]  80.76434 -76.01350

centroid$PC2
# [1] -0.4481258  0.4217655
# 




pca.results1 <- dplyr::left_join(pca.results, centroid, by = "Group", 
                                 suffix = c("",".cen"))

cols <- c('#37629C','#F18625')



pca_eig1 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[1,2],2)
pca_eig2 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[2,2],2)



PCA <- ggplot(pca.results1,aes(x=PC1,y=PC2))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(aes(color=Group),size = 3)+
  geom_segment(aes(xend=PC1.cen,yend=PC2.cen,color=Group),
               show.legend = F)+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', size = 1, fill = 'transparent'), 
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.5,"cm"),
        legend.box.spacing = unit(0,"cm"))+
  labs(x =  paste('PCA_1:', pca_eig1, '%'), y = paste('PCA_2:', pca_eig2, '%'), color = '', title = 'GSE98566')+#将 PCA 轴贡献度添加到坐标轴标题中
  scale_colour_manual(values = cols)+
  geom_label(data = centroid, 
             aes(label = Group, fill = Group), size = 4, 
             show.legend = FALSE,
             color="white")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "top")


ggpreview(PCA,width = 3,height = 3.5)
ggsave(PCA,filename = 'PCA_GSE98566.pdf',width = 3,height = 3.5,dpi = 300)
ggsave(PCA,filename = 'PCA_GSE98566.tiff',width = 3,height = 3.5,dpi = 300)
















##############################################################
################IS data sets######################################
##############################################################


#GSE16561 data pre-processing----



##1.Extract the expression matrix----

library(GEOquery)
Eset <- getGEO("GSE16561",destdir = '//GSE16561_Rawdata/')
Data <- getGEOSuppFiles(GEO = 'GSE16561',baseDir = '//GSE16561_Rawdata/')
expr_counts <- Eset$GSE16561_series_matrix.txt.gz@assayData$exprs
Pd <- Eset$GSE16561_series_matrix.txt.gz@phenoData@data

colnames(Pd)
# [1] "title"                   "geo_accession"           "status"                  "submission_date"        
# [5] "last_update_date"        "type"                    "channel_count"           "source_name_ch1"        
# [9] "organism_ch1"            "characteristics_ch1"     "characteristics_ch1.1"   "characteristics_ch1.2"  
# [13] "molecule_ch1"            "extract_protocol_ch1"    "label_ch1"               "label_protocol_ch1"     
# [17] "taxid_ch1"               "hyb_protocol"            "scan_protocol"           "description"            
# [21] "data_processing"         "platform_id"             "contact_name"            "contact_email"          
# [25] "contact_phone"           "contact_fax"             "contact_department"      "contact_institute"      
# [29] "contact_address"         "contact_city"            "contact_state"           "contact_zip/postal_code"
# [33] "contact_country"         "supplementary_file"      "data_row_count"          "age:ch1"                
# [37] "gender:ch1"              "race:ch1"     

#View sample distribution
table(str_split(string = Pd$title,pattern = '_',simplify = T)[,2])
# Control  Stroke 
# 24      39 

library(limma)


expr_counts <- 2^expr_counts 


expr_counts[is.na(expr_counts)] <- 0

expr_counts <- normalizeBetweenArrays(expr_counts)

boxplot(expr_counts[,1:10])


#探针转化


#载入注释文件
gpl<- getGEO("GPL6883", AnnotGPL=T, destdir="F:/GSE16561_Rawdata/")
ids <- Table(gpl)[,c(1,3)]

#注释
exp <- expr_counts
exp <- as.data.frame(exp)
exp$ID <- rownames(exp)
exp <- merge(x = exp,y = ids,by = 'ID')

table(duplicated(exp$`Gene symbol`))
# FALSE  TRUE 
# 17495  6929

#随机去除重复的基因
exp_unique <- exp[!duplicated(exp$`Gene symbol`),]
dim(exp_unique)
# [1] 17495    65


#删除Gene symbol列为空值的整行
exp_unique <- exp_unique[!exp_unique$`Gene symbol` == "",]
rownames(exp_unique) <- exp_unique$`Gene symbol`
exp_unique <- exp_unique[,-c(1,ncol(exp_unique))]

expr_counts <- exp_unique

expr_counts[1:4,1:4]
# GSM416528 GSM416529 GSM416530 GSM416531
# EEF1A1  1.1648112 1.3711924 1.3659349 1.5914033
# SLC35E2 0.8981245 1.0168218 1.0810391 0.9738006
# RPS28   0.7791403 0.2305806 0.7448444 0.8277220
# IPO13   1.1160861 0.9562615 0.9602976 0.9207567



##2.Organize clinical information----


colnames(Pd)
# [1] "title"                   "geo_accession"           "status"                  "submission_date"        
# [5] "last_update_date"        "type"                    "channel_count"           "source_name_ch1"        
# [9] "organism_ch1"            "characteristics_ch1"     "characteristics_ch1.1"   "characteristics_ch1.2"  
# [13] "molecule_ch1"            "extract_protocol_ch1"    "label_ch1"               "label_protocol_ch1"     
# [17] "taxid_ch1"               "hyb_protocol"            "scan_protocol"           "description"            
# [21] "data_processing"         "platform_id"             "contact_name"            "contact_email"          
# [25] "contact_phone"           "contact_fax"             "contact_department"      "contact_institute"      
# [29] "contact_address"         "contact_city"            "contact_state"           "contact_zip/postal_code"
# [33] "contact_country"         "supplementary_file"      "data_row_count"          "age:ch1"                
# [37] "gender:ch1"              "race:ch1"     


Pd <- Pd %>% 
  mutate(Group = str_split(string = title,pattern = '_',simplify = T)[,2] %>% factor(levels = c('Control','Stroke'))) %>% 
  dplyr::select(c("age:ch1","gender:ch1","race:ch1","Group")) %>% 
  dplyr::rename("Age" = "age:ch1",
                "Gender" = "gender:ch1",
                "Race" = "race:ch1")

Pd[1:5,]
# Age Gender      Race  Group
# GSM416528  48   Male Caucasian Stroke
# GSM416529  57   Male Caucasian Stroke
# GSM416530  62   Male Caucasian Stroke
# GSM416531  68   Male Caucasian Stroke
# GSM416532  60   Male Caucasian Stroke


Exp_GSE16561_Array <- expr_counts

Pd_GSE16561 <- Pd





##3.PCA----


###**PCA**

#Set up grouping
Group <- Pd_GSE16561$Group %>% factor(levels = c('Control','Stroke'))

#Start by filtering the gene to ensure that it is expressed in at least half of the sample

table(rowSums(Exp_GSE16561_Array>=1)>=round(ncol(Exp_GSE16561_Array)/2))
# FALSE  TRUE 
# 10126  7368  

Exp_GSE16561_Array_filter <- Exp_GSE16561_Array[rowSums(Exp_GSE16561_Array>=1)>=round(ncol(Exp_GSE16561_Array)/2),]
dim(Exp_GSE16561_Array_filter)
# [1] 7368   63



#Principal component analysis
set.seed(10086)
dat = as.matrix(t(Exp_GSE16561_Array_filter))

dat = dat[,colSums(is.na(dat))<nrow(dat)]
pca<-prcomp(dat,scale. = T)
pca.results <- data.frame(pca$x)[,1:2]
pca.results$Group<- Group
pca.results[1:4,1:3]
# PC1         PC2  Group
# GSM416528 -10.29247  -0.8296562 Stroke
# GSM416529 -68.27548  -0.9285190 Stroke
# GSM416530 -71.49350   0.5760536 Stroke
# GSM416531 -54.47969 -15.0376590 Stroke


centroid <- aggregate(cbind(PC1,PC2) ~ Group,
                      data = pca.results,
                      FUN = mean)

levels(centroid$Group)
# [1] "Control" "Stroke"

centroid$PC1
# [1]  4.306726 -2.650293

centroid$PC2
# [1]  22.18684 -13.65344
# 
# #因此这里我们每个组只有三个样本，因此我们对于样本的中心点需要调整，避免中心点遮挡样本点
# centroid$PC1 <- c(-50, 50)
# centroid$PC2 <- c(-5, 10)



pca.results1 <- dplyr::left_join(pca.results, centroid, by = "Group", 
                                 suffix = c("",".cen"))

cols <- c('#37629C','#F18625')



pca_eig1 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[1,2],2)
pca_eig2 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[2,2],2)



PCA <- ggplot(pca.results1,aes(x=PC1,y=PC2))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(aes(color=Group),size = 3)+
  geom_segment(aes(xend=PC1.cen,yend=PC2.cen,color=Group),
               show.legend = F)+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', size = 1, fill = 'transparent'), 
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.5,"cm"),
        legend.box.spacing = unit(0,"cm"))+
  labs(x =  paste('PCA_1:', pca_eig1, '%'), y = paste('PCA_2:', pca_eig2, '%'), color = '', title = 'GSE16561')+#将 PCA 轴贡献度添加到坐标轴标题中
  scale_colour_manual(values = cols)+
  geom_label(data = centroid, 
             aes(label = Group, fill = Group), size = 4, 
             show.legend = FALSE,
             color="white")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "top")


ggpreview(PCA,width = 3,height = 3.5)
ggsave(PCA,filename = here('/','PCA_GSE16561_All.pdf'),width = 3,height = 3.5)
ggsave(PCA,filename = here('/','PCA_GSE16561_All.tiff'),width = 3,height = 3.5,dpi = 300)

#GSE22255 data pre-processing----

library(vroom)
library(tidyverse)

##1.Extract the expression matrix----

library(GEOquery)
Eset <- getGEO("GSE22255",destdir = 'C:/GSE22255')
expr_counts <- Eset$GSE22255_series_matrix.txt.gz@assayData$exprs
Pd <- Eset$GSE22255_series_matrix.txt.gz@phenoData@data

gse_1 <- getGEO(filename = "GPL570.annot.gz")
gse <- getGEO(filename = "GPL570.soft.gz")
colnames(Pd)
# [1] "title"                     "geo_accession"             "status"                   
# [4] "submission_date"           "last_update_date"          "type"                     
# [7] "channel_count"             "source_name_ch1"           "organism_ch1"             
# [10] "characteristics_ch1"       "characteristics_ch1.1"     "characteristics_ch1.2"    
# [13] "characteristics_ch1.3"     "characteristics_ch1.4"     "characteristics_ch1.5"    
# [16] "characteristics_ch1.6"     "characteristics_ch1.7"     "characteristics_ch1.8"    
# [19] "characteristics_ch1.9"     "molecule_ch1"              "extract_protocol_ch1"     
# [22] "label_ch1"                 "label_protocol_ch1"        "taxid_ch1"                
# [25] "hyb_protocol"              "scan_protocol"             "description"              
# [28] "data_processing"           "platform_id"               "contact_name"             
# [31] "contact_email"             "contact_department"        "contact_institute"        
# [34] "contact_address"           "contact_city"              "contact_state"            
# [37] "contact_zip/postal_code"   "contact_country"           "supplementary_file"       
# [40] "data_row_count"            "age:ch1"                   "bdi:ch1"                  
# [43] "bdins:ch1"                 "bmi:ch1"                   "comorbidity:ch1"          
# [46] "education (years):ch1"     "gender:ch1"                "history of depression:ch1"
# [49] "SD:ch1"              "race:ch1"     

#View sample distribution
table(Pd$`affected status (disease state):ch1`)
#control IS patient 
#20         20 


library(limma)

expr_counts <- normalizeBetweenArrays(expr_counts)

boxplot(expr_counts[,1:10])


#No probe conversion is required





##2.Organize clinical information----


colnames(Pd)
# [1] "title"                     "geo_accession"             "status"                   
# [4] "submission_date"           "last_update_date"          "type"                     
# [7] "channel_count"             "source_name_ch1"           "organism_ch1"             
# [10] "characteristics_ch1"       "characteristics_ch1.1"     "characteristics_ch1.2"    
# [13] "characteristics_ch1.3"     "characteristics_ch1.4"     "characteristics_ch1.5"    
# [16] "characteristics_ch1.6"     "characteristics_ch1.7"     "characteristics_ch1.8"    
# [19] "characteristics_ch1.9"     "molecule_ch1"              "extract_protocol_ch1"     
# [22] "label_ch1"                 "label_protocol_ch1"        "taxid_ch1"                
# [25] "hyb_protocol"              "scan_protocol"             "description"              
# [28] "data_processing"           "platform_id"               "contact_name"             
# [31] "contact_email"             "contact_department"        "contact_institute"        
# [34] "contact_address"           "contact_city"              "contact_state"            
# [37] "contact_zip/postal_code"   "contact_country"           "supplementary_file"       
# [40] "data_row_count"            "age:ch1"                   "bdi:ch1"                  
# [43] "bdins:ch1"                 "bmi:ch1"                   "comorbidity:ch1"          
# [46] "education (years):ch1"     "gender:ch1"                "history of depression:ch1"
# [49] "SD:ch1"              "race:ch1"     


Pd <- Pd %>% 
  as.data.frame() %>%
  mutate(Group = ifelse(`insomnia:ch1` == 'control', 'Control', 'IS') %>% factor(levels = c('Control','IS'))) %>% 
  dplyr::select(c("age:ch1","bdi:ch1","bmi:ch1","comorbidity:ch1",          
                  "education (years):ch1","gender:ch1","history of depression:ch1",
                  "race:ch1","Group")) %>% 
  dplyr::rename("Age" = "age:ch1",
                "BDI" = "bdi:ch1",
                "BMI" = "bmi:ch1",
                "Comorbidity" = "comorbidity:ch1",
                "Education" = "education (years):ch1",
                "Gender" = "gender:ch1",
                "Depression" = "history of depression:ch1",
                "Race" = "race:ch1")
Pd <- Pd %>% 
  as.data.frame() %>%mutate(Group = ifelse(`affected status (disease state):ch1` == 'control', 'Control', 'IS') %>% factor(levels = c('Control','IS')))
# Age BDI         BMI Comorbidity Education Gender Depression      Race Group
# GSM6360934  65  13 21.49923325 0.638977647        16 female        yes     white    SD
# GSM6360935  75   7 26.41070366  0.95846647        16   male         no     white    SD
# GSM6360936  77   4 31.28330994  1.91693294        15 female        yes     white    SD
# GSM6360937  64   7  25.7443676           0        16 female         no non-white    SD
# GSM6360938  60   0 31.59882355           0        16   male        yes     white    SD

Exp_GSE22255_Array <- expr_counts

Pd_GSE22255 <- Pd



##3.PCA----


###**PCA**

#Set up grouping
Group <- Pd_GSE22255$Group %>% factor(levels = c('Control','IS'))

#Start by filtering the gene to ensure that it is expressed in at least half of the sample

table(rowSums(Exp_GSE22255_Array>=1)>=round(ncol(Exp_GSE22255_Array)/2))
#TRUE 
#54675   

Exp_GSE22255_Array_filter <- Exp_GSE22255_Array[rowSums(Exp_GSE22255_Array>=1)>=round(ncol(Exp_GSE22255_Array)/2),]
dim(Exp_GSE22255_Array_filter)
#[1] 54675    40


#Principal component analysis
set.seed(10086)
dat = as.matrix(t(Exp_GSE22255_Array_filter))

dat = dat[,colSums(is.na(dat))<nrow(dat)]
pca<-prcomp(dat,scale. = T)
pca.results <- data.frame(pca$x)[,1:2]
pca.results$Group<- Group
pca.results[1:4,1:3]
# PC1       PC2   Group
#GSM554014 -138.81579  19.51078 Control
#GSM554015 -135.44231  36.70737 Control
#GSM554016  -99.05315 -25.84992 Control
#GSM554017  -39.41167 -64.10688 Control


centroid <- aggregate(cbind(PC1,PC2) ~ Group,
                      data = pca.results,
                      FUN = mean)

levels(centroid$Group)
# [1] "Control"  "IS"

centroid$PC1
# [1]  20.1393 -20.1393

centroid$PC2




pca.results1 <- dplyr::left_join(pca.results, centroid, by = "Group", 
                                 suffix = c("",".cen"))

cols <- c('#37629C','#F18625')

library(FactoMineR)

pca_eig1 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[1,2],2)
pca_eig2 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[2,2],2)



PCA <- ggplot(pca.results1,aes(x=PC1,y=PC2))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(aes(color=Group),size = 3)+
  geom_segment(aes(xend=PC1.cen,yend=PC2.cen,color=Group),
               show.legend = F)+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', size = 1, fill = 'transparent'), 
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.5,"cm"),
        legend.box.spacing = unit(0,"cm"))+
  labs(x =  paste('PCA_1:', pca_eig1, '%'), y = paste('PCA_2:', pca_eig2, '%'), color = '', title = 'GSE22255')+#将 PCA 轴贡献度添加到坐标轴标题中
  scale_colour_manual(values = cols)+
  geom_label(data = centroid, 
             aes(label = Group, fill = Group), size = 4, 
             show.legend = FALSE,
             color="white")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "top")
PCA

ggpreview(PCA,width = 3,height = 3.5)
ggsave(PCA,filename = 'PCA_GSE22255_All.pdf',width = 3,height = 3.5,dpi = 300)
ggsave(PCA,filename = here('/','PCA_GSE22255_All.tiff'),width = 3,height = 3.5,dpi = 300)





#According to the PCA data, some samples were deleted
Exp_GSE22255_Array_filter <- Exp_GSE22255_Array_filter[,!colnames(Exp_GSE22255_Array_filter)%in%rownames(pca.results[pca.results$PC1>0 & pca.results$Group=='Control',])] %>% 
  as.data.frame() %>% 
  na.omit()

dim(Exp_GSE22255_Array_filter)
# [1] 33209    33

Pd_GSE22255_filter <- Pd_GSE22255[match(colnames(Exp_GSE22255_Array_filter),rownames(Pd_GSE22255)),]




set.seed(10086)
dat = as.matrix(t(Exp_GSE22255_Array_filter))

dat = dat[,colSums(is.na(dat))<nrow(dat)]
pca<-prcomp(dat,scale. = T)
pca.results <- data.frame(pca$x)[,1:2]
pca.results$Group<- Pd_GSE22255_filter$Group
pca.results[1:4,1:3]
# PC1       PC2 Group
# GSM6360934 -72.36380   3.66478    SD
# GSM6360935 -70.86557 104.14033    SD
# GSM6360936 -65.94930 110.43308    SD
# GSM6360937 -85.86780  77.00190    SD


centroid <- aggregate(cbind(PC1,PC2) ~ Group,
                      data = pca.results,
                      FUN = mean)

levels(centroid$Group)
# [1] "Control"  "SD"

centroid$PC1
# [1]  80.76434 -76.01350

centroid$PC2
# [1] -0.4481258  0.4217655





pca.results1 <- dplyr::left_join(pca.results, centroid, by = "Group", 
                                 suffix = c("",".cen"))

cols <- c('#37629C','#F18625')



pca_eig1 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[1,2],2)
pca_eig2 <- round(PCA(as.data.frame(dat),graph = FALSE)$eig[2,2],2)



PCA <- ggplot(pca.results1,aes(x=PC1,y=PC2))+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(aes(color=Group),size = 3)+
  geom_segment(aes(xend=PC1.cen,yend=PC2.cen,color=Group),
               show.legend = F)+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', size = 1, fill = 'transparent'), 
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(0.5,"cm"),
        legend.box.spacing = unit(0,"cm"))+
  labs(x =  paste('PCA_1:', pca_eig1, '%'), y = paste('PCA_2:', pca_eig2, '%'), color = '', title = 'GSE22255')+#将 PCA 轴贡献度添加到坐标轴标题中
  scale_colour_manual(values = cols)+
  geom_label(data = centroid, 
             aes(label = Group, fill = Group), size = 4, 
             show.legend = FALSE,
             color="white")+
  scale_fill_manual(values = cols)+
  theme(legend.position = "top")


ggpreview(PCA,width = 3,height = 3.5)
ggsave(PCA,filename = 'PCA_GSE22255.pdf',width = 3,height = 3.5,dpi = 300)
ggsave(PCA,filename = 'PCA_GSE22255.tiff',width = 3,height = 3.5,dpi = 300)


