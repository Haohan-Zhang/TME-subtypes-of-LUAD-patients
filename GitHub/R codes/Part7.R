#Part7
######F6B#####
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

######F6C#####
library(survival)
library(survminer)
library(tidyverse)
library(eoffice)

fit_F2 <- coxph(Surv(OS, Event=="Dead") ~ class+Stage+Age , data =clinic.class)
summary(fit_F2)

p <- ggforest(model = fit_F2,#
         data = clinic.class,
         main = 'Hazard ratio',  #
         cpositions = c(0.05, 0.15, 0.35),  #
         fontsize = 1, #
         noDigits = 3) #
pdf(file="F6C.pdf", onefile = FALSE)
print(p)
dev.off()

#####F6D#####
XCELL.MATRIX.GSE50081 <- read.csv("GSE50081.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
XCELL.MATRIX.GSE13213 <- read.csv("FGSE13213.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
XCELL.MATRIX.GSE11969 <- read.csv("GSE11969.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
XCELL.MATRIX.GSE42127 <- read.csv("GSE42127.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)
XCELL.MATRIX.GSE37745 <- read.csv("GSE37745.LUAD.XCELL.MATRIX.TCGA.50celltype.csv",
                                  row.names = 1,header = T,sep = ",",check.names = F,stringsAsFactors = F)

XCELL.MATRIX.GSE50081$CELLTYPE <- rownames(XCELL.MATRIX.GSE50081)
XCELL.MATRIX.GSE13213$CELLTYPE <- rownames(XCELL.MATRIX.GSE13213)
XCELL.MATRIX.GSE11969$CELLTYPE <- rownames(XCELL.MATRIX.GSE11969)
XCELL.MATRIX.GSE37745$CELLTYPE <- rownames(XCELL.MATRIX.GSE37745)
XCELL.MATRIX.GSE42127$CELLTYPE <- rownames(XCELL.MATRIX.GSE42127)

range(XCELL.MATRIX.GSE50081$GSM1213670)
range(XCELL.MATRIX.GSE13213$GSM333673)
range(XCELL.MATRIX.GSE11969$GSM302998)
range(XCELL.MATRIX.GSE37745$GSM1019139)
range(XCELL.MATRIX.GSE42127$GSM1032881)

XCELL.MATRIX <- merge(XCELL.MATRIX.GSE50081,XCELL.MATRIX.GSE13213,by.x="CELLTYPE",by.y="CELLTYPE")
XCELL.MATRIX <- merge(XCELL.MATRIX,XCELL.MATRIX.GSE11969,by.x="CELLTYPE",by.y="CELLTYPE")
XCELL.MATRIX <- merge(XCELL.MATRIX,XCELL.MATRIX.GSE37745,by.x="CELLTYPE",by.y="CELLTYPE")
XCELL.MATRIX <- merge(XCELL.MATRIX,XCELL.MATRIX.GSE42127,by.x="CELLTYPE",by.y="CELLTYPE")
rownames(XCELL.MATRIX) <- XCELL.MATRIX$CELLTYPE
XCELL.MATRIX <- XCELL.MATRIX[,-1]


class <-  read.csv('GSE42127+GSE37745+GSE13213+GSE11969+GSE50081.csv',
                   header = T,sep = ",",check.names = F)
head(class)

rownames(class) <- class$ID
#clinic.class$age.class <- cut(clinic.class$Age, c(-Inf , 50 , 60 , 70 , Inf) , label=c("40","50","60","70"))
#clinic.class$Smoke <- cut(clinic.class$Smoking_history, c(-Inf , 1 , Inf) , label=c("NO","YES"))
anntation.col = as.data.frame(class[,c(8,4:5,2)])
colnames(anntation.col) <- c("Class","AJCC.Stage","Gender","Age")
class(anntation.col$Class)
anntation.col$Class <- as.factor(anntation.col$Class)
anntation.col$Class <- as.numeric(anntation.col$Class)
anntation.col <- anntation.col[order(anntation.col[,1]),] #????һ?е???????
anntation.col$Class <- as.factor(anntation.col$Class)

data <- as.data.frame(t(XCELL.MATRIX))
data <- data[which(rownames(data) %in% class$ID),]

library(car)
z_norm <- function(x){
  -1 + 2*(x-min(x))/(max(x)-min(x))
}
data<-apply(data,2,z_norm)

heat.data <- as.data.frame(data)
heat.data$ID <- rownames(heat.data)
heat.data$ID

heat.data <- merge(heat.data, class,  by.x = 'ID',by.y= 'ID')
heat.data <- heat.data[order(heat.data[,58]),] 
rownames(heat.data) <- heat.data$ID
heat.data <- as.data.frame(heat.data[,c(2:51)])
heat.data <- t(heat.data)

x <- rownames(heat.data)
heat.data <- as.data.frame(apply(heat.data,2,as.numeric))
rownames(heat.data) <- x

library(pheatmap)
p <- pheatmap(heat.data, annotation_col = anntation.col,
              cluster_cols = F,
              cluster_rows = F,
              show_colnames = F,
              gaps_col = c(100,223,453),
              scale = "none",
              color = c(colorRampPalette(colors = c("#3C5488","#F8F8FF"))(40),colorRampPalette(colors = c("#F8F8FF","#E64B35"))(40)),
              #clustering_distance_rows = "manhattan", 
              # color = c(colorRampPalette(colors = c("#6998C6","#FEFEC0"))(40),colorRampPalette(colors = c("#FEFEC0","#D73027"))(40)),
              #legend_breaks=seq(-8,8,2),
              #breaks=bk
              clustering_method = "ward.D"
) 

pdf(file="F6D.pdf",width = 10, height = 8, onefile = FALSE)
print(p)
dev.off()
