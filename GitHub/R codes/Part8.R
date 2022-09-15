#Part 8 
#####F7A#####

library(pRRophetic)
library(ggplot2)
library(cowplot)
library(pheatmap)

GDSC.TCGA <- read.csv("TCGAdrugic50.icgc.csv",sep=",", header = T,row.names = 1, stringsAsFactors = F, check.names = F)
head(GDSC.TCGA)
range(GDSC.TCGA[,2:5])
GDSC.GEO <- read.csv("GEOdrug.csv",sep=",", header = T, stringsAsFactors = F, check.names = F)
head(GDSC.GEO)
range(GDSC.GEO[,2:5])


rownames(GDSC.TCGA) <- GDSC.TCGA$Drug
rownames(GDSC.GEO) <- GDSC.GEO$Drug

GDSC.TCGA <- GDSC.TCGA[,-1]
GDSC.GEO <- GDSC.GEO[,-1]

head(GDSC.TCGA)
head(GDSC.GEO)

z_norm <- function(x){
  -1 + 2*(x-min(x))/(max(x)-min(x))
}

x.TCGA <- apply(GDSC.TCGA,1,z_norm) #1行，2列
x.GEO  <- apply(GDSC.GEO,1,z_norm)

x.TCGA <- as.data.frame(t(x.TCGA))
x.GEO <- as.data.frame(t(x.GEO))

x.TCGA$DRUG <- rownames(x.TCGA)
x.GEO$DRUG <- rownames(x.GEO)
head(x.TCGA)

gdsc.heat <- merge(x.TCGA,x.GEO,by.x = "DRUG",by.y = "DRUG")

rownames(gdsc.heat) <- gdsc.heat$DRUG
head(gdsc.heat)
gdsc.heat <- gdsc.heat[,-1]
gdsc.heat$C2.x
D <- gdsc.heat[,c("C1.x","C1.y",
                  "C2.x","C2.y",
                  "C3.x","C3.y",
                  "C4.x","C4.y")] # 

annotation.col <- read.csv("TMEcluster.csv",sep=",", header = T,row.names = 1, stringsAsFactors = F, check.names = F)

p <- pheatmap(D, 
              #annotation_row =annotation.row,
              annotation_col =annotation.col,
              cluster_cols =F,
              cluster_rows = T,
              show_colnames = T,show_rownames = T,
              gaps_col = c(2,4,6),
              scale = "none",
              #clustering_distance_rows = "manhattan", 
              color = c(colorRampPalette(colors = c("#3C5488","#F8F8FF"))(40),colorRampPalette(colors = c("#F8F8FF","#E64B35"))(40)),
              #legend_breaks=seq(-8,8,2),
              #breaks=bk
              clustering_method = "ward.D",
              #annotation_colors = ann_colors,
              fontsize_col=6,
              fontsize_row=6,
              #labels_row = annotation.row$cancertype
              main ="predictive sensitivity to chemotherapy or targettherapy in GDSC "
)

p

pdf(file="7A.pdf",width = 5,height = 7, onefile = FALSE)
print(p)
dev.off()


#####F7B#####
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggthemes)
library(ggsci)

TIDE<- read.table("TIDE_value.txt",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
head(TIDE)
class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                  header = F, row.names = 1,check.names = F)
head(class)
class$ID  <- rownames(class)
colnames(class) <- c("class","ID")

TIDE <- as.data.frame(TIDE)
TIDE$ID  <- rownames(TIDE)
TIDE$ID <- substr(TIDE$ID,1,12)
#TIDE <- t(TIDE)
#colnames(TIDE) <- substr(colnames(TIDE),1,12)
#colnames(TIDE)
#TIDE$ID  <- rownames(TIDE)
clinic.TIDE <- merge(TIDE, class,  by.x = 'ID',by.y= 'ID')
head(clinic.TIDE)
#TIDE$ID <- rownames(TIDE)
class( clinic.TIDE$class)
clinic.TIDE$class <- as.factor(clinic.TIDE$class)

kruskal.test(clinic.TIDE$TIDE ~ clinic.TIDE$class)


input3e <- clinic.TIDE[,c(2,3)]
input3e$class<-as.factor(input3e$class)
input3e$TIDE <- as.numeric(input3e$TIDE)
class(input3e$TIDE)
p = 
  ggplot(input3e,aes(x=class, y=TIDE, fill=class, palette = "Lancet"), 
         ylab="TIDE score",
         xlab="", ggtheme = theme_few(),legeng="top") + 
  geom_violin( trim = F,draw_quantiles = T)+
  geom_boxplot(width = 0.2)+rotate_x_text(60)+theme_few()

p1 <- p+stat_compare_means(aes(group=class),
                           comparisons = list( c('1','2'),c('1','3'),c('1','4'),
                                               c('2','3'),c('2','4'),
                                               c('3','4')), 
                           symnum.args=list(cutpoints =c("***"=0.001, "**"=0.01, "*"=0.05)))

pdf(file="TIDE.pdf",width = 5, height = 4, onefile = FALSE)
print(p1)
dev.off()

write.csv(input3e,file="TCGA.TIDE.csv")

#####7C#####
IPS<- read.table("IPS.txt",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
head(IPS)
IPS$ID  <- rownames(IPS)
IPS$ID <- substr(IPS$ID,1,12)
clinic.IPS <- merge(IPS, class,  by.x = 'ID',by.y= 'ID')
head(clinic.IPS)
input3e <- clinic.IPS[,c(7,8)]
input3e$class<-as.factor(input3e$class)
str(input3e)

p = 
  ggplot(input3e,aes(x=class, y=IPS,fill=class,palette = "Lancet"), 
         scale_fill_brewer(""),
         ylab="IPS score",
         xlab="", ggtheme = theme_few(),legeng="top") + 
  geom_violin( trim = F,draw_quantiles = T)+
  geom_boxplot(width = 0.2)+rotate_x_text(60)+theme_few()

p1 <- p+stat_compare_means(aes(group=class),
                           comparisons = list( c('1','2'),c('1','3'),c('1','4'),
                                               c('2','3'),c('2','4'),
                                               c('3','4')),
                           symnum.args=list(cutpoints =c("***"=0.001, "**"=0.01, "*"=0.05)))

pdf(file="IPS.pdf",width = 5, height = 4, onefile = FALSE)
print(p1)
dev.off()

write.csv(input3e,file="TCGA.IPS.csv")
