#Part 2

###################FIG1###################
#####FIG1A#####
library(ggplot2)
setwd("F:\\肺癌第一次返修\\F1")
XCELL.MATRIX <- read.csv("FXCELL.MATRIX.TCGA+GSE68465.csv",
                         header = T, row.names = 1,check.names = F)
XCELL.MATRIX[1:6,1:6]
clinic <- read.csv("GSE68465.TCGA.热图.clinic.csv",
                   header = T, row.names = 1,check.names = F)
clinic[1:6,1:6]

x <- as.data.frame(t(XCELL.MATRIX))
x$ID <- rownames(x)
clinic$ID <- rownames(clinic)
clinical.data <- merge(x,clinic,by.x="ID",by.y = "ID")
write.csv(clinical.data,"clinical.data.csv")


clinical.data <- read.csv("clinical.data.csv",row.names = 1,header = T,sep=",", check.names = F)
#spearman
corResult_test=as.data.frame(apply(clinical.data[,2:51],2,function(x){
  cor.test(x,as.numeric(clinical.data$Age) ,method="spearman")$estimate
}))
corResult_test.p=as.data.frame(apply(clinical.data[,2:51],2,function(x){
  cor.test(x,clinical.data$Age,method="spearman")$p.value
}))
corResult <- cbind(corResult_test, corResult_test.p)
colnames(corResult) <- c("cor", "p")
corResult$celltype<-rownames(corResult)


p<- ggplot()+ geom_point(data = corResult, 
                         mapping = aes(x=celltype,y=1,size = -log10(p),color = cor),
                         alpha=0.85, 
                         group = "Class")+
  scale_color_gradient2(low = "#5D90BA",
                        mid = "white",
                        high = "#D20A13",
                        breaks=seq(-1,1,0.2))+
  theme_bw() +  
  theme_classic()+ 
  theme(
    # Rotate the x-axis lables so they are legible
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),#??ȥy????ֵ
    axis.ticks.y = element_blank(),#??ȥy???̶???
    # Force the plot into a square aspect ratio
    aspect.ratio = 0.3)

pdf(file="1A.pdf",width = 18,height = 5, onefile = FALSE)
print(p)
dev.off()

corResult_test$celltype <- rownames(corResult_test)
corResult_test.p$celltype <- rownames(corResult_test.p)

colnames(corResult_test) <- c("correlation","celltype")
colnames(corResult_test.p) <- c("P","celltype")

y <- merge(corResult_test,corResult_test.p,by.x="celltype",by.y = "celltype")
write.csv(y,"age.csv")


#####FIG1B#####
clinical.data <- read.csv('clinical.data.csv',sep=",", header = T,row.names = 1, stringsAsFactors = F, check.names = F)
z_norm<-function(x){
  (x-min(x))/(max(x)-min(x))
}

outTab2=data.frame()

clinical.data$sex <- ifelse(clinical.data$Sex  =="Female",0,1)
clinical.data$tnm <- ifelse(clinical.data$AJCC_Stage =="Stage I/II",0,1)
clinical.data$t <- ifelse(clinical.data$T_Stage =="T1/T2",0,1)
clinical.data$n <- ifelse(clinical.data$N =="N0",0,1)
clinical.data$m <- ifelse(clinical.data$M =="M0",0,1)

for (i in colnames(clinical.data[,c(2:51)])) {
  clinical.data$z.norm <- z_norm(clinical.data[,i])
  Fc <- mean(clinical.data[which(clinical.data$m == 1),i])/
    mean(clinical.data[which(clinical.data$m == 0),i])
  p <- t.test(clinical.data[,i] ~ clinical.data$m )$p.value
  outTab2 = rbind(outTab2, cbind(i,Fc, p))
} 
outTab2$Fc <- as.character(outTab2$Fc)
outTab2$Fc <- as.numeric(outTab2$Fc)
outTab2$p <- as.character(outTab2$p)
outTab2$p <- as.numeric(outTab2$p)
outTab2$reg <- ifelse(outTab2$Fc >1 , 'Up', 'Down')
outTab2$text <- ifelse((outTab2$Fc >1 ) & outTab2$p < 0.05, 1, 
                       ifelse((outTab2$Fc < 1 ) & outTab2$p < 0.05, 2, 0))
#outTab2$text <- factor(outTab2$text,levels=c(0,1,2))
p.m<- 
  ggplot(data = outTab2, aes(x= -log10(p), y = Fc )) +
  guides(colour=guide_legend(title=NULL))+
  #theme(legend.position = "none")+
  geom_point(data = outTab2,aes(x= -log10(p), y = Fc, color= factor(text)),size = 3)+
  scale_color_manual(values = mycol)+ #ylim(0,3)+
  theme_bw()+ geom_hline( yintercept = 1) + geom_vline(xintercept = -log10(0.05))+
  ggrepel::geom_text_repel(data = subset(outTab2, outTab2$p < 0.05),size= 4,aes(label= i ),
                           box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"), 
                           segment.color = "black", show.legend = FALSE)

library(cowplot)
library(ggplot2)
plot_grid(p.t)

plot_grid(p.sex, p.tnm,p.t,p.n,p.m, nrow = 5)

pdf(file="1Bp.m.pdf",width = 6,height = 6, onefile = FALSE)
print(p.m)
dev.off()

#####FIG1C#####
library(reshape2)
library(corrplot)
library(plyr)
library(igraph) 
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

poscol <- "#FB9A99" #
negcol <- "#C6DBEF" #
mycol <- c("#FDBF6F", "#1F78B4", "#E31A1C", "#8C510A") 

XCELL.MATRIX <- read.csv("XCELL.MATRIX.TCGA+GSE68465.csv",
                         row.names =1, sep = "," , header = T,check.names = F, stringsAsFactors = T)
XCELL.MATRIX[1:6,1:6]
clinic<-read.csv('clinical.data.csv',header = T,sep = ",",check.names = F)
clinic[1:6,1:6]

x<-as.data.frame(t(XCELL.MATRIX))
input_data <- x[which(rownames(x) %in% clinic$ID),]

data <- XCELL.MATRIX[which(rownames(XCELL.MATRIX) %in% clinic$ID),]


clinical.data
gg <- data.frame(c(1:50),1,1,1,1,1,1,1)
colnames(gg) <- c("cell_types","coef","HR","se","z","p","ll","ul")
for(i in 2:51){
  ll <- coxph(Surv(OS, Event == "Dead") ~ clinical.data[,i], data =  clinical.data)
  beta <- coef(ll)
  se <- sqrt(diag(vcov(ll)))
  HR <- exp(beta)
  HRse <- HR * se
  z = beta/se
  p = 1 - pchisq((beta/se)^2, 1)
  HRCILL = exp(beta - qnorm(.975, 0, 1) * se)
  HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)
  gg[i-1,1] <- colnames(clinical.data)[i]
  gg[i-1,2] <- beta
  gg[i-1,3] <- HR
  gg[i-1,4] <- se
  gg[i-1,5] <- z
  gg[i-1,6] <- p
  gg[i-1,7] <- HRCILL
  gg[i-1,8] <- HRCIUL
  
}
write.csv(gg,file = "50celltype_cox_LUAD.csv")

##
bb <- read.csv("50celltype_cox_LUAD.csv",header = T,row.names = 1);
head(bb)
bb$cell_types <- as.character(bb$cell_types) 
rownames(bb) <- bb$cell_types
head(bb)
#bb <- bb[,-1]
colnames(bb)[1] <- c("ID")
#??pvalue???ƽڵ?Բ?Ĵ?С
bb$weight <- abs(log10(bb$p))
#??HR??Բ?ĵ?????ɫ
bb$weight_HR <- (as.numeric(bb$HR)-1)*100
bb$colr <- ifelse(bb$weight_HR<0, "green", "black")
head(bb)

head(input_data)
head(corr)
#单因素生存分析HR
write.csv(bb, "output_HR_corlorx.csv", quote = F)


corr <- cor(input_data, method = "spearman")


cor.mtest <- function(corr, ...) {
  corr <- as.matrix(corr)
  n <- ncol(corr)
  p.corr <- matrix(NA, n, n)
  diag(p.corr) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(corr[, i],method = "spearman", corr[, j], ...)
      p.corr[i, j] <- p.corr[j, i] <- tmp$p.value
    }
  }
  colnames(p.corr) <- rownames(p.corr) <- colnames(corr)
  p.corr
}

p.corr <- cor.mtest(input_data) 
head(p.corr[, 1:5])

#?ϲ?????ϵ????Pֵ
rr <- as.data.frame(corr);
rr$ID <- rownames(rr)
cor <- melt(rr,"ID",value.name = "cor"); #head(cor)

pp <- as.data.frame(p.corr);
pp$ID <- rownames(pp)
pvalue <- melt(pp,"ID",value.name = "pvalue"); #head(pvalue)
colnames(pvalue) <- c("from","to","pvalue")

corpvlue <- cbind(pvalue, cor)
head(corpvlue)
corpvlue <- corpvlue[, -c(4:5)]
head(corpvlue)
dim(corpvlue)
corpvlue$fdr <- p.adjust(corpvlue$pvalue)


corpvlue <- corpvlue[!corpvlue$cor==1,]
dim(corpvlue)

corpvlue$weight <- corpvlue$pvalue
corpvlue$weight <- -log10(corpvlue$weight)
summary(duplicated(corpvlue$weight))
corpvlue <- corpvlue[!duplicated(corpvlue$weight),] 
dim(corpvlue)
#write.csv(corpvlue,"corpvlue.csv",sep=",")

corpvlue <- corpvlue[which(corpvlue$fdr < 0.0001 ),] #
corpvlue <- corpvlue[-which(abs(corpvlue$cor) <0.5 ),] #
dim(corpvlue)

head(corpvlue)

corpvlue$color <- ifelse(corpvlue$cor<0, negcol, poscol)



nodes <- read.csv("nodes.csv",header = T)
head(nodes)
summary(nodes$media %in% bb$ID) #
nodes <- merge(nodes, bb, by.x = "media", "ID", all.x = T, all.y = T) #

nodes$Fraction <- abs(nodes$weight_HR)

nodes$id <- paste("S", 01:50, sep = "") #
nodes <- nodes[order(nodes$type),]
nodes <- nodes[,c(ncol(nodes),1:ncol(nodes)-1)] 
nodes <- nodes[order(nodes$type),]
head(nodes)


celltoS <- as.vector(nodes[,c(1)])
names(celltoS) <- nodes[,c(2)]
#names(celltoS) <- paste0(nodes[,c(2)], "_XCELL")
corpvlue$from <- revalue(as.character(corpvlue$from) ,celltoS)
corpvlue$to <- revalue(as.character(corpvlue$to) ,celltoS)
(links <- corpvlue)

#
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
# Generate colors based on cell clusters:
V(net)$color <- revalue(nodes$type,c("Lymphoids"=mycol[1],
                                     "Myeloids"=mycol[2],
                                     "Stem Cells"=mycol[3],
                                     "Stromal cells"=mycol[4]))
# Compute node degrees (#links) and use that to set node size:
# Set edge width based on weight-log10(p_value):
V(net)$size <- (1 + V(net)$weight)*3.8 #

V(net)$label <- V(net)$media #
E(net)$arrow.mode <- 0 #
E(net)$edge.color <- "tomato" # tomato gray80
E(net)$width <- 1+E(net)$weight/15  

plot(net,
     layout=layout_in_circle, #
     edge.curved=.2, #
     vertex.label.color=V(net)$color, 
     vertex.label.dist= -2, 
     edge.color=links$color)

legend("topright",#
       inset=-0.3,
       c("Lymphoids", "Myeloids", "Stem Cells","Stromal cells"),
       pch=21, 
       col="black", pt.bg=mycol, pt.cex=3,
       cex=1.3, bty="n", ncol=1)


f <- c(0.1, 0.05, 0.01, 0.001)
s <- (1 + -log10(f))
legend("bottomright", 
       inset=c(0,-.25), 
       legend=f, text.width = .2, 
       title = "logrank test, P value", title.adj = -.5,
       pch=21, pt.cex=s, bty='n',
       horiz = TRUE, #
       col = "black")


legend("bottomright",
       c("Positive correlation with FDR < 0.0001", 
         "Negative correlation with FDR < 0.0001"),
       inset=c(0,-.2),
       col = c(poscol, negcol), bty="n", 
       cex = 1, lty = 1, lwd = 5)


#####Fig1D#####
library(survival)
library(ggplot2)
library(forestplot)
library(ggplot2)

clinical.data <- read.csv('clinical.data.csv',sep=",", header = T,row.names = 1, stringsAsFactors = F, check.names = F)
z_norm<-function(x){
  (x-min(x))/(max(x)-min(x))
}
outTab=data.frame()

class(clinical.data$Age)
clinical.data$Age <- as.numeric(clinical.data$Age)
# + age + as.factor(AJCC.Stage)  + A18_Sex;+ Race + Smoking_history + A18_Sex
for (i in colnames(clinical.data[,c(2:51)])) {
  fit<-coxph(Surv(OS, Event=="Dead")~ z_norm(clinical.data[,i]) 
             +as.numeric(Age) + as.factor(AJCC.Stage) 
             ,data=clinical.data,x=T,y=T)
  coxSummary = summary(fit)
  id=row.names(coxSummary$coefficients)
  coef=coxSummary$coefficients[,"coef"]
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  HR=coxSummary$conf.int[,"exp(coef)"]
  HR.95L=coxSummary$conf.int[,"lower .95"]
  HR.95H=coxSummary$conf.int[,"upper .95"]
  Tab=cbind(Cell.types=i,coef,HR, HR.95L,HR.95H,coxP)
  outTab=rbind(outTab,Tab)
} 
write.csv(outTab,"50cell.survival.multi.withRC.csv")


forest.uni <- read.csv("50cell.survival.multi.withRC.csv",header = T);

par(mar=(c(10,10,9,8)))
mycol <- rep(c("#000000", "#E69F00", "#0072B2", "#009E73", "#D55E00",
               "#223D6C","#D20A13","#FFD121","#088247","#11AA4D",
               "#58CDD9","#7A142C","#5D90BA",
               "#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D",
               "#7CC767" ,"#999999","#FF0099",  "#56B4E9",
               "#F0E442", "#CC79A7","#990000","#9900cc",
               "#66FF66","#663300","#0000FF","#CC0033","#FF0000","#000099",
               "#660066","#333333","#FFCCCC","#993366","#33CC33","#CC9900" ),4)

class(forest.uni$HR)
x<-forest.uni$Cell.types
forest.uni <- as.data.frame(apply(forest.uni[,2:6],2,as.numeric))
forest.uni$Cell.types <- x

forest.uni<-forest.uni[order(forest.uni$HR,decreasing = T),]

p <- forestplot(as.matrix(forest.uni[,1]),#   
                forest.uni$HR, forest.uni$HR.95L , forest.uni$HR.95H ,
                graph.pos= 2, graphwidth = unit(60,"mm"), 
                #is.summary = F,fn.ci_norm ="fpDrawBarCI",#
                lwd.ci = 2,lty.ci = 1,
                col=fpColors(box=c("#BC3C29FF","steelblue","#0072B5FF"), lines=c("black",'black'), 
                             zero = "black"),#
                ci.vertices.height = 0.1, # 
                boxsize=0.3 ,
                clip = c(0,0),
                # legend = 'Subgroup analysis for risk stratification',
                txt_gp=fpTxtGp(ticks = gpar(cex=1),summary = gpar(cex=1),cex=1),
                line.margin=unit(10,"mm"),
                lineheight = unit(4,'mm'),  
                vertices = TRUE,xlab="Hazard Ratio",
                #hrzl_lines=list("2" = gpar(lwd=2, col="black"), #
                #               "4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#
                #             "38" = gpar(lwd=2, col="black")),#
                colgap = unit(3,"mm"), zero=1,xticks = c(0,0.25,0.5,1.0,2,3,5))

pdf(file="1D.pdf",width = 9,height = 18, onefile = FALSE)
print(p)
dev.off()

