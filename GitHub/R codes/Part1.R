library(immunedeconv)
library(xCell)
######celltype######
type <- c(
  #Lymphoids
  "CD4+ memory T-cells"  ,        
  "CD4+ naive T-cells" ,          
  "CD4+ T-cells"   ,              
  "CD4+ Tcm"  ,                    
  "CD4+ Tem"  ,                   
  "CD8+ naive T-cells"  ,         
  "CD8+ T-cells"  ,               
  "CD8+ Tcm"  ,                    
  "CD8+ Tem"   ,  
  "B-cells"    , 
  "Tgd cells"     ,    
  "Th1 cells" ,
  "Th2 cells" , 
  "Tregs",
  "Memory B-cells"   , 
  "Class-switched memory B-cells" ,
  #"naive B-cells",
  "NK cells"    ,                  
  "NKT" ,
  "pro B-cells"   ,
  "Plasma cells" ,         
  #Stem Cells
  #"MEP",                    
  "GMP"     ,   #Granulocyte-macrophage progenitor
  "HSC",        #Hematopoietic stem cells 
  "MPP"     ,   #Multipotent progenitors
  "CLP"     ,   #common lymphoid progenitor   
  "CMP" ,       #common myeloid progenito
  #"Erythrocytes",           
  #"Platelets"   ,        
  #Myeloids
  "Basophils"    ,          
  "cDC"     ,    #Conventional dendritic cells
  "DC"    ,      #dendritic cell 
  "Eosinophils"    ,        
  "iDC"    ,     #Immature dendritic cells
  "Macrophages"            ,      
  "Macrophages M1"   ,             
  "Macrophages M2"  ,             
  "Mast cells"   , 
  "pDC" ,        #Plasmacytoid dendritic cells
  "Neutrophils" ,            
  "Monocytes"  ,            
  "aDC" ,
  #Stromal cells
  "Chondrocytes"      ,     
  "Endothelial cells"     , 
  "Fibroblasts",           
  "ly Endothelial cells"  , 
  "mv Endothelial cells" ,  
  "Megakaryocytes"  ,       
  #"Smooth muscle",         
  "Pericytes"  ,               
  "Preadipocytes" ,         
  "Adipocytes"   ,         
  "MSC"  ,       #Mesenchymal stem cells 
  "Osteoblast" ,            
  #"Skeletal muscle"  ,     
  "Myocytes"               
  #Others
  #"Astrocytes"  ,            
  #"Epithelial cells"    ,   
  #"Keratinocytes"  ,        
  #"Sebocytes",             
  #"Neurons"   ,               
  #"Mesangial cells"   ,      
  #"Hepatocytes" ,           
  #"Melanocytes"          
)
#####survival#####
library(sva)
library(ConsensusClusterPlus)
library(survminer)
library(survival)
library(estimate)
library(pheatmap)
library(ggthemes)
library(tableone)
library(rms) 
library(ggplot2)


expr.tpm.tcga <- read.table("uniq.tumor.tpm.txt",
                            row.names = 1,header = T,check.names = F,stringsAsFactors = F)
expr.tpm.tcga[1:6,1:6]

range(expr.tpm.tcga) 
expr.tpm.tcga <- log2(expr.tpm.tcga + 1)
range(expr.tpm.tcga) # 0-20

expr.tpm.GSE68465 <- read.csv("GSE68465.genematrix.csv",
                              row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
rownames(expr.tpm.GSE68465) <- expr.tpm.GSE68465$Symbol
expr.tpm.GSE68465 <- expr.tpm.GSE68465[,-1]
expr.tpm.GSE68465[1:6,1:6]

range(expr.tpm.GSE68465) 
expr.tpm.GSE68465 <- log2(expr.tpm.GSE68465 + 1)
range(expr.tpm.GSE68465) # 0-20
expr.set.GSE68465 <- expr.tpm.GSE68465


comgene <- intersect(intersect(rownames(expr.tpm.tcga), rownames(expr.tpm.GSE68465)))
combined.expr <- cbind.data.frame(tcga.expr[comgene,],
                                  gse41613.expr[comgene,],
                                  gse65858.expr[comgene,])

expr.set.GSE68465 <- expr.set.GSE68465[which(rownames(expr.set.GSE68465) %in% df1$Symbol),]
expr.tpm.tcga <- expr.tpm.tcga[which(rownames(expr.tpm.tcga) %in% df1$Symbol),]
expr.set.GSE68465$gene <- rownames(expr.set.GSE68465)
expr.tpm.tcga$gene <- rownames(expr.tpm.tcga)
combined.expr <- merge(expr.tpm.tcga,expr.set.GSE68465,by.x="gene",by.y="gene")
rownames(combined.expr) <- combined.expr$gene
combined.expr[1:6,1:6]
combined.expr <- combined.expr[,-which(colnames(combined.expr) %in% c("gene"))]
expr.tpm.tcga <- expr.tpm.tcga[,-which(colnames(expr.tpm.tcga) %in% c("gene"))]
expr.set.GSE68465 <- expr.set.GSE68465[,-which(colnames(expr.set.GSE68465) %in% c("gene"))]


batch <- data.frame(batch = rep(c("TCGA","GSE68465"), times = c(ncol(expr.tpm.tcga),ncol(expr.set.GSE68465))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- ComBat(dat=as.matrix(combined.expr), 
                               batch=batch$batch, 
                               mod=modcombat)

combined.expr.combat <- read.csv("combined.expr.combat.csv",
                                 row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
clinic <- read.csv('GSE68465.TCGA.clinic.插补后.csv',header = T,sep = ",",check.names = F)


combined.expr.combat <- combined.expr.combat[,which(colnames(combined.expr.combat) %in% clinic$ID)]

XCELL.MATRIX <- xCellAnalysis(combined.expr.combat,
                              rnaseq = F,
                              cell.types.use =type,
                              parallel.type = "SOCK")
XCELL.MATRIX.PVALUE <- xCellSignifcanceBetaDist(XCELL.MATRIX, 
                                                beta_params = NULL, 
                                                rnaseq = F,
                                                file.name = NULL)


write.csv(combined.expr.combat,"combined.expr.combat.csv")
write.csv(XCELL.MATRIX,"XCELL.MATRIX.TCGA+GSE68465.csv")
write.csv(XCELL.MATRIX.PVALUE,"XCELL.MATRIX.PVALUE.TCGA+GSE68465.csv")


XCELL.MATRIX <- read.csv('XCELL.MATRIX.TCGA+GSE68465.csv',
                         header = T,sep = ",",check.names = F)
data <- as.data.frame(t(XCELL.MATRIX))
data <- data[which(rownames(data) %in% clinic$ID),]

library(car)
z_norm <- function(x){
  -1 + 2*(x-min(x))/(max(x)-min(x))
}
data<-apply(data,2,z_norm)

data<-t(data)
data <- as.matrix(data)
getwd()
res2 <- ConsensusClusterPlus(data, 
                             maxK = 6, #
                             reps = 1000, 
                             pItem = 0.7, 
                             pFeature = 1,
                             title="scale50cell3",
                             clusterAlg = 'km', #k-means
                             distance="euclidean", 
                             seed=1234,
                             plot="pdf", #??"png"
                             corUse="all.obs",writeTable=T)#pairwise.complete.obs
class.anno <- as.data.frame(res2[[4]]$consensusClass,row.names = names(res2[[4]]$consensusClass) ) 
class.anno$ID <- rownames(class.anno)
colnames(class.anno) <- c('class', "ID")
class(class.anno$class)
class.anno$class <- as.numeric(class.anno$class)
table(class.anno)

clinic <- read.csv('GSE68465.TCGA.clinic.csv',header = T,sep = ",",check.names = F)
data <- as.data.frame(t(XCELL.MATRIX))
data <- data[which(rownames(data) %in% clinic$ID),]

clinic.class <- merge(clinic, class.anno , by.x = "ID", by.y = "ID")
clinic.class$class
clinic.class$class<-as.factor(clinic.class$class)
clinic.class$OS<-as.numeric(clinic.class$OS)

clinic.class$class  <- ifelse(clinic.class$class  == 1, "B",
                              ifelse(clinic.class$class  ==2, "A",
                                     ifelse(clinic.class$class  == 3, "C","D"
                                     )))  

fit <- coxph(Surv(OS, Event=="Dead") ~ class , data =clinic.class)
summary(fit)

head(clinic.class)
class(clinic.class$Age)
clinic.class$Age <- as.numeric(clinic.class$Age)

fit2 <- coxph(Surv(OS, Event=="Dead") ~ class+AJCC.Stage , data =clinic.class)
summary(fit2)
fit3 <- coxph(Surv(OS, Event=="Dead") ~ class+Age+AJCC.Stage , data =clinic.class)
summary(fit3)

fit1 <- survfit(Surv(OS, Event=="Dead") ~ class , data =clinic.class)

p <- ggsurvplot(
  fit1,   data = clinic.class,   
  pval = TRUE, pval.size = 3,   #
  ggtheme = theme_few(),       
  conf.int = F,  #
  legend.title = "TMErisk",  #
  xlim = c(0,120*30),  ylim = NULL, #
  xlab = "Time in years",  
  break.time.by = 24*30,  # 
  xscale=12*30,
  risk.table =T, # 
  surv.median.line = "hv", #
  risk.table.y.text.col = T,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 5,
  # legend.labs =  c("high","Low" ),   
  palette = "nejm", # 
  font.y  = c(12, "bold"), 
  font.x  = c(12, "bold"),
  legend = 'top',  #
  font.legend = c(8, "bold"))

p
pdf(file="2D.pdf",width = 8,height = 8, onefile = FALSE)
print(p)
dev.off()

######neuralnet######
library(caret)
library(neuralnet)
library(dplyr)
library(pROC)

class.anno <- read.csv("scale50cell3.k=4.consensusClass.csv",header = F)
colnames(class.anno) <- c("ID","class")
head(class.anno)

XCELL.MATRIX <- read.csv("XCELL.MATRIX.TCGA+GSE68465.csv", 
                         row.names = 1,header = T,sep = ",",check.names = F)
cell <- as.data.frame(t(XCELL.MATRIX))
cell$ID <- rownames(cell)
data <- merge(cell,class.anno,by="ID")
rownames(data) <- data$ID
mydata <- data[,-1]

set.seed(1234)

class(mydata$class)
mydata$class <- as.factor(mydata$class)
index <- createDataPartition(mydata$class, p=0.7, list=F)
traindata <- mydata[index, ]
testdata <- mydata[-index, ]


standard <- preProcess(traindata, method = c("center","scale"))
traindata_std <- predict(standard, traindata)
#standard <- preProcess(testdata, method = c("center","scale"))
testdata_std <- predict(standard, testdata)

traindatax <- traindata_std[,-51]
traindatax$one = traindata_std$class == 1
traindatax$two = traindata_std$class == 2
traindatax$three = traindata_std$class == 3
traindatax$four = traindata_std$class == 4

colnames(traindatax)
coln <- colnames(traindatax)
coln <- gsub("[-]","",coln)
coln <- gsub("[ ]","",coln)
coln <- gsub("[+]","",coln)
colnames(traindatax) <- coln
n <- names(traindatax)
n

f <- as.formula(paste("one + two + three + four ~", paste(n[!n %in% c("one" ,"two","three","four")], collapse = " + "))) 
bp_modal <- neuralnet(f ,
                      data=traindatax,
                      threshold=0.01,
                      stepmax=100000,
                      rep = 1,
                      err.fct="sse",
                      linear.output= F,
                      hidden= c(12))

testdata_std_x <- testdata_std[,-51]
coln <- colnames(testdata_std_x)
coln
coln <- gsub("[-]","",coln)
coln <- gsub("[ ]","",coln)
coln <- gsub("[+]","",coln)
colnames(testdata_std_x) <- coln
colnames(testdata_std_x)

bp_result <- predict(bp_modal,testdata_std_x)

#predict.table
#predict.table = table(testdata_std$class,bp_pred)
#predict.table

bp_pred =c("one","two","three","four")[apply(bp_result,1, which.max)]%>%as.data.frame()
bp_pred_num <- bp_pred %>% mutate(pred=case_when(
  bp_pred=='one' ~ 1,
  bp_pred=='two' ~ 2,
  bp_pred=='three' ~ 3,
  bp_pred=='four' ~ 4
))

bp_pred1 <- bp_pred_num[[2]] %>% as.factor()

bp_roc <- multiclass.roc(testdata_std$class,as.numeric(bp_pred1),direction = "<")
bp_roc

getwd()
save(bp_modal, file = "bp_modal0.9585.RData")


library(survival)
library(survminer)
library(tableone)
library(ggthemes)
library(ggplot2)
library(rms) 
load("bp_modal0.9585.RData")
plot(bp_modal)

####GSE135222#####
expr.set.GSE135222 <- read.csv("GSE135222.genematrix.csv",
                               row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
expr.set.GSE135222 <- expr.set.GSE135222[,-1]
clinic <- read.csv('GSE135222.clinic.csv',
                   header = T,sep = ",",check.names = F)
expr.set.GSE135222 <- expr.set.GSE135222[which(colnames(expr.set.GSE135222) %in% clinic$ID)]
XCELL.MATRIX <- xCellAnalysis(expr.set.GSE135222,
                              rnaseq = T,
                              cell.types.use =type,
                              parallel.type = "SOCK")
XCELL.MATRIX.PVALUE <- xCellSignifcanceBetaDist(XCELL.MATRIX, 
                                                beta_params = NULL, 
                                                rnaseq = T,
                                                file.name = NULL)
getwd()
write.csv(XCELL.MATRIX,"GSE135222.LUAD.XCELL.MATRIX.TCGA.50celltype.csv")
write.csv(XCELL.MATRIX.PVALUE,"GSE135222.LUAD.XCELL.MATRIX.PVALUE.50celltype.csv")

XCELL.MATRIX <- read.csv("GSE135222.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                        row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
XCELL.MATRIX.PVALUE <- read.csv("GSE135222.LUAD.XCELL.MATRIX.PVALUE.50celltype.csv",
                               row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)



x<- as.matrix(t(XCELL.MATRIX))
standard <- preProcess(x, method = c("center","scale"))
GEO_testdata_std <- predict(standard, x)

colnames(GEO_testdata_std)
coln <- colnames(GEO_testdata_std)
coln <- gsub("[ ]","",coln)
coln <- gsub("[+]","",coln)
coln <- gsub("[-]","",coln)
colnames(GEO_testdata_std) <- coln

bp_result <- predict(bp_modal,GEO_testdata_std)
#bp_result$net.result

bp_pred =c("one","two","three","four")[apply(bp_result,1, which.max)]%>%as.data.frame()
bp_pred_num <- bp_pred%>%mutate(pred=case_when(
  bp_pred=='one' ~ 1,
  bp_pred=='two' ~ 2,
  bp_pred=='three' ~ 3,
  bp_pred=='four' ~ 4
))

bp_pred1 <- bp_pred_num[[2]] %>% as.factor()
rownames(bp_pred_num) <- rownames(GEO_testdata_std)
bp_pred_num

validation.class<- as.data.frame(bp_pred_num)
validation.class$ID <- rownames(validation.class)
validation.class <- validation.class[,-1]
colnames(validation.class) <- c("class","ID")
head(validation.class)

clinic.geo.class <- merge(clinic, validation.class , by.x = "ID", by.y = "ID")
clinic.geo.class$class

class(clinic.geo.class$class)
clinic.geo.class$class<-as.factor(clinic.geo.class$class)

head(clinic.geo.class)

clinic.geo.class$class  <- ifelse(clinic.geo.class$class  == 4, "A",
                                  ifelse(clinic.geo.class$class  == 1, "B",
                                         ifelse(clinic.geo.class$class  == 2, "C","D"
                                         )))

fit <- coxph(Surv(`pfs.time`, `progression-free survival (pfs)`=='1') ~ class, data =clinic.geo.class)
summary(fit)

class(clinic.geo.class$age)
fit2 <- coxph(Surv(`pfs.time`, `progression-free survival (pfs)`=='1') ~ class + age, data =clinic.geo.class)
summary(fit2)

fit1 <- survfit(Surv(`pfs.time`, `progression-free survival (pfs)`=='1') ~ class , data =clinic.geo.class)

p <- ggsurvplot(
  fit1,   data = clinic.geo.class,   
  pval = TRUE, pval.size = 3,   
  ggtheme = theme_few(),       
  conf.int = F,  
  legend.title = "TMErisk",  
  xlim = c(0,24*30),  ylim = NULL, 
  xlab = "Time in months",  
  break.time.by = 12*30,  
  xscale=12*30,
  risk.table =T, 
  surv.median.line = "hv", 
  risk.table.y.text.col = T,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 5,
  # legend.labs =  c("high","Low" ),   
  palette = "nejm", 
  font.y  = c(12, "bold"), 
  font.x  = c(12, "bold"),
  legend = 'top', 
  font.legend = c(8, "bold"))

p
pdf(file="7D.pdf",width = 8,height = 8, onefile = FALSE)
print(p)
dev.off()

write.csv(clinic.geo.class,"GSE135222.clinic.class.csv")

######GSE42127+GSE37745+GSE13213+GSE11969+GSE50081######
library(survival)
library(survminer)
library(tableone)
library(ggthemes)
library(ggplot2)
library(rms) 
load("bp_modal0.9079.RData")
plot(bp_modal)

#50081
XCELL.MATRIX.GSE50081 <- read.csv("GSE50081.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
clinic <- read.csv('GSE50081.LUAD.clinic.csv',
                   header = T,sep = ",",check.names = F)
XCELL.MATRIX.GSE50081 <- XCELL.MATRIX.GSE50081[which(colnames(XCELL.MATRIX.GSE50081) %in% clinic$geo_accession)]
class(XCELL.MATRIX.GSE50081)

x<- as.matrix(t(XCELL.MATRIX.GSE50081))
standard <- preProcess(x, method = c("center","scale"))
GEO_testdata_std <- predict(standard, x)

colnames(GEO_testdata_std)
coln <- colnames(GEO_testdata_std)
coln <- gsub("[ ]","",coln)
coln <- gsub("[+]","",coln)
coln <- gsub("[-]","",coln)
colnames(GEO_testdata_std) <- coln

bp_result <- predict(bp_modal,GEO_testdata_std)
#bp_result$net.result

#????Ԥ??ֵ
bp_pred =c("one","two","three","four")[apply(bp_result,1, which.max)]%>%as.data.frame()
bp_pred_num <- bp_pred%>%mutate(pred=case_when(
  bp_pred=='one' ~ 1,
  bp_pred=='two' ~ 2,
  bp_pred=='three' ~ 3,
  bp_pred=='four' ~ 4
))

bp_pred1 <- bp_pred_num[[2]] %>% as.factor()
rownames(bp_pred_num) <- rownames(GEO_testdata_std)
bp_pred_num

validation.class<- as.data.frame(bp_pred_num)
validation.class$ID <- rownames(validation.class)
validation.class <- validation.class[,-1]
colnames(validation.class) <- c("class","ID")
head(validation.class)

clinic$ID <- clinic$geo_accession
clinic.geo.class <- merge(clinic, validation.class , by.x = "ID", by.y = "ID")
clinic.geo.class$class

class(clinic.geo.class$class)
clinic.geo.class$class<-as.factor(clinic.geo.class$class)

head(clinic.geo.class)

fit <- coxph(Surv(`survival time:ch1`, `status:ch1`=='dead') ~ class, data =clinic.geo.class)
summary(fit)
write.csv(clinic.geo.class,"GSE50081.clinic.class.csv")

#GSE42127
XCELL.MATRIX.GSE42127 <- read.csv("GSE42127.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
clinic <- read.csv('\GSE42127.LUAD.clinic.csv',
                   header = T,sep = ",",check.names = F)
XCELL.MATRIX.GSE42127 <- XCELL.MATRIX.GSE42127[which(colnames(XCELL.MATRIX.GSE42127) %in% clinic$geo_accession)]
class(XCELL.MATRIX.GSE42127)

x<- as.matrix(t(XCELL.MATRIX.GSE42127))
standard <- preProcess(x, method = c("center","scale"))
GEO_testdata_std <- predict(standard, x)

colnames(GEO_testdata_std)
coln <- colnames(GEO_testdata_std)
coln <- gsub("[ ]","",coln)
coln <- gsub("[+]","",coln)
coln <- gsub("[-]","",coln)
colnames(GEO_testdata_std) <- coln

bp_result <- predict(bp_modal,GEO_testdata_std)
#bp_result$net.result

bp_pred =c("one","two","three","four")[apply(bp_result,1, which.max)]%>%as.data.frame()
bp_pred_num <- bp_pred%>%mutate(pred=case_when(
  bp_pred=='one' ~ 1,
  bp_pred=='two' ~ 2,
  bp_pred=='three' ~ 3,
  bp_pred=='four' ~ 4
))

bp_pred1 <- bp_pred_num[[2]] %>% as.factor()
rownames(bp_pred_num) <- rownames(GEO_testdata_std)
bp_pred_num

validation.class<- as.data.frame(bp_pred_num)
validation.class$ID <- rownames(validation.class)
validation.class <- validation.class[,-1]
colnames(validation.class) <- c("class","ID")
head(validation.class)

clinic$ID <- clinic$geo_accession
clinic.geo.class <- merge(clinic, validation.class , by.x = "ID", by.y = "ID")
clinic.geo.class$class

class(clinic.geo.class$class)
clinic.geo.class$class<-as.factor(clinic.geo.class$class)

head(clinic.geo.class)

fit <- coxph(Surv(`overall_survival_months`, `Event`=='D') ~ class, data =clinic.geo.class)
summary(fit)
write.csv(clinic.geo.class,"GSE42127.clinic.class.csv")

#GSE37745
XCELL.MATRIX.GSE37745 <- read.csv("GSE37745.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
clinic <- read.csv('GSE37745.LUAD.clinic.csv',
                   header = T,sep = ",",check.names = F)
XCELL.MATRIX.GSE37745 <- XCELL.MATRIX.GSE37745[which(colnames(XCELL.MATRIX.GSE37745) %in% clinic$geo_accession)]
class(XCELL.MATRIX.GSE37745)

x<- as.matrix(t(XCELL.MATRIX.GSE37745))
standard <- preProcess(x, method = c("center","scale"))
GEO_testdata_std <- predict(standard, x)

colnames(GEO_testdata_std)
coln <- colnames(GEO_testdata_std)
coln <- gsub("[ ]","",coln)
coln <- gsub("[+]","",coln)
coln <- gsub("[-]","",coln)
colnames(GEO_testdata_std) <- coln

bp_result <- predict(bp_modal,GEO_testdata_std)
#bp_result$net.result

bp_pred =c("one","two","three","four")[apply(bp_result,1, which.max)]%>%as.data.frame()
bp_pred_num <- bp_pred%>%mutate(pred=case_when(
  bp_pred=='one' ~ 1,
  bp_pred=='two' ~ 2,
  bp_pred=='three' ~ 3,
  bp_pred=='four' ~ 4
))

bp_pred1 <- bp_pred_num[[2]] %>% as.factor()
rownames(bp_pred_num) <- rownames(GEO_testdata_std)
bp_pred_num

validation.class<- as.data.frame(bp_pred_num)
validation.class$ID <- rownames(validation.class)
validation.class <- validation.class[,-1]
colnames(validation.class) <- c("class","ID")
head(validation.class)

clinic$ID <- clinic$geo_accession
clinic.geo.class <- merge(clinic, validation.class , by.x = "ID", by.y = "ID")
clinic.geo.class$class

class(clinic.geo.class$class)
clinic.geo.class$class<-as.factor(clinic.geo.class$class)

head(clinic.geo.class)

clinic.geo.class$class  <- ifelse(clinic.geo.class$class  == 4, "A",
                                  ifelse(clinic.geo.class$class  == 1, "B",
                                         ifelse(clinic.geo.class$class  == 2, "C","D"
                                         )))

fit <- coxph(Surv(`days to determined death status:ch1`, `dead:ch1`=='yes') ~ class, data =clinic.geo.class)
summary(fit)

write.csv(clinic.geo.class,"GSE37745.clinic.class.csv")

#GSE13213
XCELL.MATRIX.GSE13213 <- read.csv("GSE13213.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
clinic <- read.csv('GSE13213.clinic.csv',
                   header = T,sep = ",",check.names = F)
XCELL.MATRIX.GSE13213 <- XCELL.MATRIX.GSE13213[which(colnames(XCELL.MATRIX.GSE13213) %in% clinic$geo_accession)]
class(XCELL.MATRIX.GSE13213)

x<- as.matrix(t(XCELL.MATRIX.GSE13213))
standard <- preProcess(x, method = c("center","scale"))
GEO_testdata_std <- predict(standard, x)

colnames(GEO_testdata_std)
coln <- colnames(GEO_testdata_std)
coln <- gsub("[ ]","",coln)
coln <- gsub("[+]","",coln)
coln <- gsub("[-]","",coln)
colnames(GEO_testdata_std) <- coln

bp_result <- predict(bp_modal,GEO_testdata_std)
#bp_result$net.result

bp_pred =c("one","two","three","four")[apply(bp_result,1, which.max)]%>%as.data.frame()
bp_pred_num <- bp_pred%>%mutate(pred=case_when(
  bp_pred=='one' ~ 1,
  bp_pred=='two' ~ 2,
  bp_pred=='three' ~ 3,
  bp_pred=='four' ~ 4
))

bp_pred1 <- bp_pred_num[[2]] %>% as.factor()
rownames(bp_pred_num) <- rownames(GEO_testdata_std)
bp_pred_num

validation.class<- as.data.frame(bp_pred_num)
validation.class$ID <- rownames(validation.class)
validation.class <- validation.class[,-1]
colnames(validation.class) <- c("class","ID")
head(validation.class)

clinic$ID <- clinic$geo_accession
clinic.geo.class <- merge(clinic, validation.class , by.x = "ID", by.y = "ID")
clinic.geo.class$class

class(clinic.geo.class$class)
clinic.geo.class$class<-as.factor(clinic.geo.class$class)

head(clinic.geo.class)

fit <- coxph(Surv(`Survival (days):ch1`, `Status:ch1`=='Dead') ~ class, data =clinic.geo.class)
summary(fit)
write.csv(clinic.geo.class,"GSE13213.clinic.class.csv")

#GSE11969
XCELL.MATRIX.GSE11969 <- read.csv("GSE11969.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
clinic <- read.csv('GSE11969.LUAD.clinic.csv',
                   header = T,sep = ",",check.names = F)
XCELL.MATRIX.GSE11969 <- XCELL.MATRIX.GSE11969[which(colnames(XCELL.MATRIX.GSE11969) %in% clinic$geo_accession)]
class(XCELL.MATRIX.GSE11969)

x<- as.matrix(t(XCELL.MATRIX.GSE11969))
standard <- preProcess(x, method = c("center","scale"))
GEO_testdata_std <- predict(standard, x)

colnames(GEO_testdata_std)
coln <- colnames(GEO_testdata_std)
coln <- gsub("[ ]","",coln)
coln <- gsub("[+]","",coln)
coln <- gsub("[-]","",coln)
colnames(GEO_testdata_std) <- coln

bp_result <- predict(bp_modal,GEO_testdata_std)
#bp_result$net.result

#
bp_pred =c("one","two","three","four")[apply(bp_result,1, which.max)]%>%as.data.frame()
bp_pred_num <- bp_pred%>%mutate(pred=case_when(
  bp_pred=='one' ~ 1,
  bp_pred=='two' ~ 2,
  bp_pred=='three' ~ 3,
  bp_pred=='four' ~ 4
))

bp_pred1 <- bp_pred_num[[2]] %>% as.factor()
rownames(bp_pred_num) <- rownames(GEO_testdata_std)
bp_pred_num

validation.class<- as.data.frame(bp_pred_num)
validation.class$ID <- rownames(validation.class)
validation.class <- validation.class[,-1]
colnames(validation.class) <- c("class","ID")
head(validation.class)

clinic$ID <- clinic$geo_accession
clinic.geo.class <- merge(clinic, validation.class , by.x = "ID", by.y = "ID")
clinic.geo.class$class

class(clinic.geo.class$class)
clinic.geo.class$class<-as.factor(clinic.geo.class$class)

head(clinic.geo.class)

fit <- coxph(Surv(`Survival (days):ch1`, `Status:ch1`=='Dead') ~ class, data =clinic.geo.class)
summary(fit)
write.csv(clinic.geo.class,"GSE11969.clinic.class.csv")


clinic.class <-  read.csv('GSE42127+GSE37745+GSE13213+GSE11969+GSE50081.csv',
                          header = T,sep = ",",check.names = F)

class(clinic.class$class)
clinic.class$class<-as.factor(clinic.class$class)

head(clinic.class)

clinic.class$class  <- ifelse(clinic.class$class  == 3, "A",
                              ifelse(clinic.class$class  == 1, "B",
                                     ifelse(clinic.class$class  == 2, "C","D"
                                     )))

fit <- coxph(Surv(OS, Event=='Dead') ~ class, data =clinic.class)
summary(fit)

fit <- coxph(Surv(OS, Event=='Dead') ~ class+Stage+Age, data =clinic.class)
summary(fit)

fit1 <- survfit(Surv(OS, Event=='Dead') ~ class , data =clinic.class)

p <- ggsurvplot(
  fit1,   data = clinic.class,   
  pval = TRUE, pval.size = 3,   #
  ggtheme = theme_few(),       
  conf.int = F,  #
  legend.title = "TMErisk",  #
  xlim = c(0,120*30),  ylim = NULL, #
  xlab = "Time in years",  
  break.time.by = 24*30,  # 
  xscale=12*30,
  risk.table =T, # 
  surv.median.line = "hv", #
  risk.table.y.text.col = T,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 5,
  # legend.labs =  c("high","Low" ),   
  palette = "nejm", # 
  font.y  = c(12, "bold"), 
  font.x  = c(12, "bold"),
  legend = 'top',  #
  font.legend = c(8, "bold"))
p 


pdf(file="F6B.pdf",width = 7,height = 6, onefile = FALSE)
print(p)
dev.off()




