#Part4
######3A######
XCELL.MATRIX <- read.csv("XCELL.MATRIX.TCGA+GSE68465.csv",
                         header = T, row.names = 1,check.names = F)
clinic <- read.csv("GSE68465.TCGA.clinic.csv",
                   header = T, row.names = 1,check.names = F)
class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                  header = F, row.names = 1,check.names = F)
colnames(class) <-  ("Class")
class$ID <- rownames(class)
clinic$ID <- rownames(clinic)

clinic.class <- merge(class,clinic,by.x = "ID",by.y = "ID")

rownames(clinic.class) <- clinic.class$ID
#clinic.class$age.class <- cut(clinic.class$Age, c(-Inf , 50 , 60 , 70 , Inf) , label=c("40","50","60","70"))
#clinic.class$Smoke <- cut(clinic.class$Smoking_history, c(-Inf , 1 , Inf) , label=c("NO","YES"))
anntation.col = as.data.frame(clinic.class[,c(2,9,5:7,11,3)])
colnames(anntation.col) <- c("Class","AJCC.Stage","T.Stage","N.Stage","M.Stage","Gender","Age")
class(anntation.col$Class)
anntation.col$Class <- as.factor(anntation.col$Class)
anntation.col$Class <- as.numeric(anntation.col$Class)
anntation.col <- anntation.col[order(anntation.col[,1]),] #????һ?е???????
anntation.col$Class <- as.factor(anntation.col$Class)

data <- as.data.frame(t(XCELL.MATRIX))
data <- data[which(rownames(data) %in% clinic$ID),]

library(car)
z_norm <- function(x){
  -1 + 2*(x-min(x))/(max(x)-min(x))
}
data<-apply(data,2,z_norm)

heat.data <- as.data.frame(data)
heat.data$ID <- rownames(heat.data)
heat.data$ID

heat.data <- merge(heat.data, class,  by.x = 'ID',by.y= 'ID')
heat.data <- heat.data[order(heat.data[,52]),] 
rownames(heat.data) <- heat.data$ID
heat.data <- as.data.frame(heat.data[,c(2:51)])
heat.data <- t(heat.data)

x <- rownames(heat.data)
heat.data <- as.data.frame(apply(heat.data,2,as.numeric))
rownames(heat.data) <- x

library(pheatmap)
p <- pheatmap(heat.data, annotation_col = anntation.col,
              cluster_cols = F,
              show_colnames = F,
              gaps_col = c(157,407,752),
              scale = "none",
              color = c(colorRampPalette(colors = c("#3C5488","#F8F8FF"))(40),colorRampPalette(colors = c("#F8F8FF","#E64B35"))(40)),
              #clustering_distance_rows = "manhattan", 
              # color = c(colorRampPalette(colors = c("#6998C6","#FEFEC0"))(40),colorRampPalette(colors = c("#FEFEC0","#D73027"))(40)),
              #legend_breaks=seq(-8,8,2),
              #breaks=bk
              clustering_method = "ward.D"
) 

pdf(file="3A",width = 10, height = 8, onefile = FALSE)
print(p)
dev.off()

######3B######
#EPIC
library(EPIC)
library(dplyr)

expr.tpm <- read.csv("combined.expr.combat.csv",
                     row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
range(expr.tpm)

out <- EPIC(bulk = expr.tpm,  #a matrix (nGenes x nSamples) ,PM, RPKM or FPKM
            reference = TRef,  #"BRef" or "TRef" or referenceCellsList
            #mRNA_cell = mRNA_cell_vector, 
            #sigGenes = sigGenes_vector 
) 

out$cellFractions
out$mRNAProportions
out$fit.gof

cellFractions <- out$cellFractions

mydata <- as.data.frame(cellFractions)
head(mydata)

clinic.class <- read.csv("scale50cell3.k=4.consensusClass.csv",header = F, row.names = 1,check.names = F)
head(clinic.class)
clinic.class$ID  <- rownames(clinic.class)
colnames(clinic.class) <- c("class","ID")

mydata$ID <- rownames(mydata)
mydata <- merge(mydata,clinic.class,by.x = 'ID',by.y = 'ID')
rownames(mydata) <- mydata$ID 
head(mydata)
mydata <- mydata[,-1]
class(mydata$Bcells)

datamean=group_by(mydata, class) %>% summarize_each(funs(mean))
write.csv(datamean,"EPIC.csv")

data <- read.csv("EPICok.csv",header = T,check.names = F)
data

p <- ggplot(data,aes(x=class,y=value,fill=Celltype)) +
  geom_bar(stat="identity")+
  theme_classic() +
  #scale_fill_brewer("Tones")+ 
  labs(x='TMECluater',y = 'EPIC')+
  scale_fill_manual(values = mycolor)+
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"))
pdf(file="3BEPIC.pdf",width = 5, height = 4, onefile = FALSE)
print(p)
dev.off()

#MCPcounter
library(curl)
library(devtools)
library(MCPcounter)

expr.tpm <- read.csv("combined.expr.combat.csv",
                     row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
range(expr.tpm)

probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),stringsAsFactors=FALSE,colClasses='character')#
genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep='\t',stringsAsFactors=FALSE,header=TRUE,colClasses='character',check.names=FALSE)

#run example code
results <- MCPcounter.estimate( 
  expr.tpm, #	matrix or data.frame with features x samples，MCPcounterExampleData
  featuresType=c('HUGO_symbols'), #'affy133P2_probesets','ENSEMBL_ID','ENTREZ_ID'  #选择分析特征，有三种可供选择（探针信息、symbol、ID，减少名称转换）
  probesets = probesets, #
  genes = genes)

head(results)

mydata <- as.data.frame(t(results))
head(mydata)

clinic.class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                         header = F, row.names = 1,check.names = F)
head(clinic.class)
clinic.class$ID  <- rownames(clinic.class)
colnames(clinic.class) <- c("class","ID")

mydata$ID <- rownames(mydata)
mydata <- merge(mydata,clinic.class,by.x = 'ID',by.y = 'ID')
rownames(mydata) <- mydata$ID 
head(mydata)
mydata <- mydata[,-1]
class(mydata$`T cells`)

datamean=group_by(mydata, class) %>% summarize_each(funs(mean))
datamean <- as.data.frame(t(datamean))
write.csv(datamean,"MCPcounter.csv")

data <- read.csv("MCPcounterok.csv",
                 header = T,check.names = F)
data

p <- ggplot(data,aes(x=class,y=value,fill=Celltype)) +
  geom_bar(stat="identity",position = "fill")+
  theme_classic() +
  #scale_fill_brewer("Tones")+ 
  labs(x='TMECluater',y = 'MCPcounter')+
  scale_fill_manual(values = mycolor)+
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"))

pdf(file="3BMCPcounter.pdf",width = 5, height = 4, onefile = FALSE)
print(p)
dev.off()

#QUANTISEQ
library(quantiseqr)
library("dplyr")
library("ggplot2")
library("tidyr")
library("tibble")
library("GEOquery")
library("reshape2")
library("SummarizedExperiment")

expr.tpm <- read.csv("combined.expr.combat.csv",
                     row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
range(expr.tpm)

ti_racle <- quantiseqr::run_quantiseq(expression_data = expr.tpm,
                                      signature_matrix = "TIL10",
                                      is_arraydata = FALSE,
                                      is_tumordata = TRUE,
                                      scale_mRNA = TRUE
)

quantiplot(ti_racle)


mydata <- as.data.frame(ti_racle)
head(mydata)

clinic.class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                         header = F, row.names = 1,check.names = F)
head(clinic.class)
clinic.class$ID  <- rownames(clinic.class)
colnames(clinic.class) <- c("class","ID")

mydata <- merge(mydata,clinic.class,by.x = 'Sample',by.y = 'ID')

head(mydata)
rownames(mydata) <- mydata$Sample 
mydata <- mydata[,-1]
class(mydata$B.cells)

datamean=group_by(mydata, class) %>% summarize_each(funs(mean))
datamean <- as.data.frame(t(datamean))
write.csv(datamean,"QUANTISEQ.csv")

data <- read.csv("QUANTISEQok.csv",
                 header = T,check.names = F)

p <- ggplot(data,aes(x=class,y=value,fill=Celltype)) +
  geom_bar(stat="identity")+
  theme_classic() +
  #scale_fill_brewer("Tones")+ 
  labs(x='TMECluater',y = 'QUANTISEQ')+
  scale_fill_manual(values = mycolor)+
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"))

pdf(file="3BQUANTISEQ.pdf",width = 5, height = 4, onefile = FALSE)
print(p)
dev.off()

#CIBERSOFT
sig_matrix <- read.table("LM22.txt", 
                         sep = "\t", check.names = F, stringsAsFactors = F, header = T, row.names = 1)

mixture_file <-  read.csv("combined.expr.combat.csv",
                          row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
range(mixture_file)

source("CIBERSORT.R")
res_cibersort <- CIBERSORT(sig_matrix, 
                           mixture_file, 
                           perm=1000, 
                           QN = T)
save(res_cibersort,file = "CIBERSORT.RData")

res_cibersort <- as.data.frame(res_cibersort)
class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                  header = F, row.names = 1,check.names = F)
head(class)
colnames(class) <-  ("Class")
class$ID <- rownames(class)
clinic$ID <- rownames(clinic)

head(res_cibersort)
data <- res_cibersort[,1:22]
data$ID <- rownames(data)
clinic.class <- merge(data,class,by.x = "ID",by.y = "ID")
head(clinic.class)

rownames(clinic.class) <- clinic.class$ID
#clinic.class$age.class <- cut(clinic.class$Age, c(-Inf , 50 , 60 , 70 , Inf) , label=c("40","50","60","70"))
#clinic.class$Smoke <- cut(clinic.class$Smoking_history, c(-Inf , 1 , Inf) , label=c("NO","YES"))
anntation.col = as.data.frame(clinic.class[,c(1,24)])
class(anntation.col$Class)
anntation.col$Class <- as.factor(anntation.col$Class)
anntation.col$Class <- as.numeric(anntation.col$Class)
anntation.col <- anntation.col[order(anntation.col[,2]),] #????һ?е???????
anntation.col$Class <- as.factor(anntation.col$Class)
row <- rownames(anntation.col)
anntation.col <- as.data.frame(anntation.col[,-1])
rownames(anntation.col) <- row
colnames(anntation.col) <- c("TMECluster")
head(anntation.col)

data <- data[which(rownames(data) %in% clinic$ID),]
data <- data[,-which(colnames(data) %in% c("ID"))]

library(car)
z_norm <- function(x){
  -1 + 2*(x-min(x))/(max(x)-min(x))
}
data<-apply(data,2,z_norm)

heat.data <- as.data.frame(data)
heat.data$ID <- rownames(heat.data)
heat.data$ID

heat.data <- merge(heat.data, class,  by.x = 'ID',by.y= 'ID')
heat.data <- heat.data[order(heat.data[,24]),] 
rownames(heat.data) <- heat.data$ID
heat.data <- as.data.frame(heat.data[,c(2:23)])
heat.data <- t(heat.data)

x <- rownames(heat.data)
heat.data <- as.data.frame(apply(heat.data,2,as.numeric))
rownames(heat.data) <- x

library(pheatmap)
p <- pheatmap(heat.data, 
              annotation_col = anntation.col,
              cluster_cols = F,
              show_colnames = F,
              gaps_col = c(157,407,752),
              scale = "none",
              color = c(colorRampPalette(colors = c("#3C5488","#F8F8FF"))(40),colorRampPalette(colors = c("#F8F8FF","#E64B35"))(40)),
              #clustering_distance_rows = "manhattan", 
              # color = c(colorRampPalette(colors = c("#6998C6","#FEFEC0"))(40),colorRampPalette(colors = c("#FEFEC0","#D73027"))(40)),
              #legend_breaks=seq(-8,8,2),
              #breaks=bk
              clustering_method = "ward.D"
) 

table(anntation.col)

pdf(file="SF1CIBERSORT.pdf",width = 10, height = 8, onefile = FALSE)
print(p)
dev.off()

#####3C#####
library(ggpubr)
library(ggthemes)
estimate <- read.table("estimate.scores.txt",
                       row.names = 1,header = T,check.names = F,stringsAsFactors = F)
head(estimate)

coln <- rownames(estimate)
coln <- gsub("[.]","-",coln)
rownames(estimate) <- coln

data <- estimate[which(rownames(estimate) %in% clinic.class$ID),]
tcga.class <- clinic.class[which((clinic.class$ID) %in% rownames(estimate)),]
data$ID <- rownames(data)

data <- merge(data,tcga.class,by.x = 'ID',by.y= 'ID')
head(data)

rownames(data) <- data$ID
class(data$class)
data$class <- as.factor(data$class)
head(data)
rt <- data[,c(5,6)]

head(rt)

p = 
  ggplot(rt,aes(x=class, y=TumorPurit,fill=class,palette = "nejm"), 
         #legend.title = "Leukocyte Fraction",
         ylab="Stromal Score",
         xlab="TME Cluster", ggtheme = theme_few(),legeng="top") + 
  geom_violin( trim = F,draw_quantiles = T)+
  geom_boxplot(width = 0.2)+rotate_x_text(60)+theme_few()

p1<- p+stat_compare_means(aes(group=pre.class),
                          comparisons = list( c('1','2'),c('1','3'),c('1','4'),
                                              c('2','3'),c('2','4'),
                                              c('3','4')), 
                          symnum.args=list(cutpoints =c("***"=0.001, "**"=0.01, "*"=0.05)))
p1

pdf(file="3C.pdf",width = 5, height = 4, onefile = FALSE)
print(p1)
dev.off()

immunesubtype <- read.csv("immunelandscape.csv",header = T,sep = ",")
rownames(immunesubtype) <- immunesubtype$ID
clinic.class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                         header = F, row.names = 1,check.names = F)
head(clinic.class)
clinic.class$ID  <- rownames(clinic.class)
colnames(clinic.class) <- c("class","ID")

data <- immunesubtype[which((immunesubtype$ID) %in% clinic.class$ID),]
tcga.class <- clinic.class[which((clinic.class$ID) %in% immunesubtype$ID),]
data <- merge(data,tcga.class,by.x = 'ID',by.y= 'ID')
head(data)

rownames(data) <- data$ID
class(data$class)
data$class <- as.factor(data$class)
data[1:6,1:10]
rt <- data[,c(7,65)]
head(rt)

p = 
  ggplot(rt,aes(x=class, y=Intratumor.Heterogeneity,fill=class,palette = "nejm"), 
         #legend.title = "Leukocyte Fraction",
         ylab="Intratumor Heterogeneity",
         xlab="TME Clusters", ggtheme = theme_few(),legeng="top") + 
  geom_violin( trim = F,draw_quantiles = T)+
  geom_boxplot(width = 0.2)+rotate_x_text(60)+theme_few()

p1<- p+stat_compare_means(aes(group=pre.class),
                          comparisons = list( c('1','2'),c('1','3'),c('1','4'),
                                              c('2','3'),c('2','4'),
                                              c('3','4')), 
                          symnum.args=list(cutpoints =c("***"=0.001, "**"=0.01, "*"=0.05)))
p1

pdf(file="3C.pdf",width = 5, height = 4, onefile = FALSE)
print(p1)
dev.off()

#####3D#####
library(ggalluvial)
library(dplyr)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

clinical.data <- read.csv("TCGA.clinical.data.csv",
                          row.names = 1,header = T,sep=",", check.names = F)

sankey <- merge(sankey,stage,by.x = "ID",by.y = "ID")
sankey$AJCC.Stage
rownames(sankey) <- sankey$ID
head(sankey)
sankey <- sankey[,-1]
#mycol <- rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)

df <- sankey
head(df)
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
df <- removeRowsAllNa(df)

UCB_lodes <- to_lodes_form(df[,1:ncol(df)],
                           key = "class",
                           axes = 1:ncol(df),
                           id = "Cohort")
dim(UCB_lodes)
head(UCB_lodes)
tail(UCB_lodes)
ggplot(UCB_lodes,
       aes(x = class, stratum = stratum, alluvium = Cohort,
           fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 1/8) + #
  geom_stratum(alpha = .9,width = 1/10) + #
  geom_text(stat = "stratum", size = 3, color="black") + #

  scale_fill_manual(values = mycol) +
  
  xlab("") + ylab("") +
  theme_bw() + #
  theme(panel.grid =element_blank()) + #
  theme(panel.border = element_blank()) + #
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + #去掉坐标轴
  ggtitle("")+
  guides(fill = FALSE) 


clinic.class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                         header = F, row.names = 1,check.names = F)
head(clinic.class)
clinic.class$ID  <- rownames(clinic.class)
colnames(clinic.class) <- c("class","ID")

immunesubtype <- read.csv("immunelandscape.csv",header = T,sep = ",")
head(immunesubtype)
immunesubtype[1:6,1:6]
rownames(immunesubtype) <- immunesubtype$ID
data <- immunesubtype[which((immunesubtype$ID) %in% clinic.class$ID),]

tcga.class <- clinic.class[which((clinic.class$ID) %in% immunesubtype$ID),]

data <- merge(data,tcga.class,by.x = 'ID',by.y= 'ID')
data[1:6,1:6]
rownames(data) <- data$ID
class(data$class)
data$class <- as.factor(data$class)

mydata <- data[,c(9:14)]

mydata[is.na(mydata)] <- min(mydata,na.rm = T)*0.01
range(mydata)

library(car)
z_norm <- function(x){
  -1 + 2*(x-min(x))/(max(x)-min(x))
}
mydata <- apply(mydata,2,z_norm)

anntation.col = as.data.frame(tcga.class[,c(1)])#[,c(2:3,6:9, 12:13,17,19)]
rownames(anntation.col) <- rownames(tcga.class)
colnames(anntation.col) <- c("class")
class(anntation.col$class)
anntation.col$class <- as.factor(anntation.col$class)

head(mydata)
class(mydata)
mydata <- as.data.frame(mydata)
mydata$ID <- rownames(mydata)
head(clinic.class)
heat.data <- merge(mydata, clinic.class,  by.x = 'ID',by.y= 'ID')
head(heat.data)
class(heat.data$class)
heat.data <- heat.data[order(heat.data[,8]),] 
rownames(heat.data) <- heat.data$ID

coln <- rownames(heat.data)
heat.data <- as.data.frame(heat.data[,c(2:7)])
rownames(heat.data) <-coln
heat.data <- t(heat.data)

#heat.data[is.na(heat.data)] <- min(heat.data,na.rm = T)*0.01
library(pheatmap)
pheatmap(heat.data, annotation_col = anntation.col,
         cluster_cols = F,
         show_colnames = F,
         gaps_col = c(79,228,404),scale = "none",
         clustering_distance_rows = "manhattan", 
         #color = c(colorRampPalette(colors = c("#6998C6","#FEFEC0"))(40),colorRampPalette(colors = c("#FEFEC0","#D73027"))(40)),
         color = c(colorRampPalette(colors = c("#3C5488","#F8F8FF"))(40),colorRampPalette(colors = c("#F8F8FF","#E64B35"))(40)),
         #legend_breaks=seq(-8,8,2),
         #breaks=bk
         clustering_method = "ward.D"
) 

#####3E#####
Immunomodulators.name <- read.csv("Immunomodulators.name.csv",
                                  header = T,sep=",", check.names = F)
head(Immunomodulators.name)

clinic <- read.csv("scale50cell3.k=4.consensusClass.csv",
                   header = F, row.names = 1,check.names = F)
head(clinic)
clinic$ID  <- rownames(clinic)
colnames(clinic) <- c("class","ID")

expr <- read.csv("combined.expr.combat.csv",
                 header = T, row.names = 1,check.names = F)
expr[1:6,1:6]

Immunomodulators.expr <- expr[which(rownames(expr) %in% Immunomodulators.name$gene),which(colnames(expr) %in% clinic$ID)]

write.csv(Immunomodulators.expr,"Immunomodulators.expr.csv")

Immunomodulators.expr <- read.csv('Immunomodulators.expr.tpm.csv', header = T, sep = ",",check.names = F,row.names = 1)


Immunomodulators.expr[1:6,1:6]
range(Immunomodulators.expr)

Immunomodulators.expr$gene <- rownames(Immunomodulators.expr)
Immunomodulators.name$gene
Immunomodulators.expr <- merge(Immunomodulators.expr,Immunomodulators.name,by.x = 'gene',by.y = 'gene')
Immunomodulators.expr <- Immunomodulators.expr[order(Immunomodulators.expr$class),] 
rownames(Immunomodulators.expr) <- Immunomodulators.expr$gene
Immunomodulators.expr <- Immunomodulators.expr[-which(colnames(Immunomodulators.expr)%in%c("gene","class"))]

data <- as.data.frame(t(Immunomodulators.expr))

data$ID <- rownames(data)
data$ID
data <- merge(data,clinic,by.x = 'ID',by.y = 'ID')
data <- data[order(data$class),] 
head(data)
rownames(data) <- data$ID
data <- data[-which(colnames(data)%in%c("ID","class"))]

library(car)
z_norm <- function(x){
  -1 + 2*(x-min(x))/(max(x)-min(x))
}
data<-apply(data,2,z_norm)

data[1:6,1:6]
data <- as.data.frame(t(data))
range(data)

annol_row = Immunomodulators.name[which(Immunomodulators.name$gene%in%rownames(data)),]
row <- annol_row$gene
annol_row <- as.data.frame(annol_row[,-1]); rownames(annol_row) <- row
annol_col = clinic[order(clinic[,1]),] 
col <- annol_col$ID
annol_col <- as.data.frame(annol_col[,-2]); rownames(annol_col) <- col
rownames(data)

p<- pheatmap(data,
             annotation_col = annol_col,
             annotation_row = annol_row,
             scale = "none",
             gaps_row = c(11,14,19,22,42,47),
             cluster_cols = F, 
             cluster_rows =  F,
             show_colnames = F,
             gaps_col = c(157,407,752),
             color = c(colorRampPalette(colors = c("#3C5488","#F8F8FF"))(40),colorRampPalette(colors = c("#F8F8FF","#E64B35"))(40)),
             #color = c(colorRampPalette(colors = c("blue","white","red"))(60)),
             #legend_breaks = seq(-8,8,2),
             #breaks= c(seq(-5,0,by=0.1),seq(0,5,by=0.1)),
             clustering_method = "average",
             border_color = "black",
             border = T
) 

pdf(file="3E.pdf",width = 8, height = 8, onefile = FALSE)
print(p)
dev.off()




