#Part3
#####2E#####
library(survival)
library(survminer)
library(tidyverse)
class.anno <- read.csv('scale50cell3.k=4.consensusClass.csv',
                       header = F,sep = ",",check.names = F)
head(class.anno)
colnames(class.anno) <- c('ID', "class")

clinic <- read.csv('GSE68465.TCGA.clinic.csv',header = T,sep = ",",check.names = F)


clinic.class <- merge(clinic, class.anno , by.x = "ID", by.y = "ID")
clinic.class$class
clinic.class$class<-as.factor(clinic.class$class)
clinic.class$OS<-as.numeric(clinic.class$OS)
fit_F2 <- coxph(Surv(OS, Event=="Dead") ~ class+Age+AJCC.Stage+Sex , data =clinic.class)
summary(fit_F2)

ggforest(model = fit_F2,#
         data = clinic.class,
         main = 'Hazard ratio',  #
         cpositions = c(0.05, 0.15, 0.35),  #
         fontsize = 1, #
         noDigits = 3) #

#####2F#####
library(limma)
library(data.table)
library(tidyr)
library(SeuratObject) 
library(monocle)
counts <- read.csv('combined.expr.combat.csv',
                   header = T,row.names = 1,sep = ",",check.names = F)
clinic.class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                         header = F,row.names = 1)
clinic.class$ID  <- rownames(clinic.class)
colnames(clinic.class) <- c("class","ID")
head(clinic.class)


class <- tidyr::unite(clinic.class, "class_ID", class, ID)
head(class)


x <- as.data.frame(t(counts))
x$ID <- rownames(x)
#head(x)
colnames(class) <- c("row")
class$ID <- rownames(class)
x <- merge(class,x,by.x = "ID",by.y = "ID")
rownames(x) <- x$row
x = x[,-which(colnames(x) %in% c("ID","row"))]

count <- as.data.frame(t(x))  
head(count)

write.table(count,"scRNA.count.txt")
save(count,file = "scRNA.count.RData")


pbmc <- CreateSeuratObject(counts = count, 
                           project = "seurat",
                           names.delim = "_")

data <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix') #

pd <- new('AnnotatedDataFrame', data = pbmc@meta.data) #
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData) #


mycds  <- newCellDataSet(data, #
                         phenoData = pd, #
                         featureData = fd, #
                         lowerDetectionLimit = 0.5, 
                         expressionFamily = negbinomial.size()) #expression response variables

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)


dfgene1<-read.table("LUAD_limma_test_result.1_vs_Others.txt",sep="\t",header=T,row.names=1)
dfgene2<-read.table("LUAD_limma_test_result.2_vs_Others.txt",sep="\t",header=T,row.names=1)
dfgene3<-read.table("LUAD_limma_test_result.3_vs_Others.txt",sep="\t",header=T,row.names=1)
dfgene4<-read.table("LUAD_limma_test_result.4_vs_Others.txt",sep="\t",header=T,row.names=1)
dfgene<-rbind(dfgene1,dfgene2,dfgene3,dfgene4)

diff.genes <- unique(rownames(subset(dfgene,FDR<0.05 & abs(log2FC)>1)))

#disp_table <- dispersionTable(mycds)
#disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds1 <- setOrderingFilter(mycds, diff.genes)
#p3 <- plot_ordering_genes(mycds)

#降维
mycds2 <- reduceDimension(mycds1, 
                          max_components = 4,  #
                          method = 'tSNE' ,  #
                          norm_method="log", #
                          pseudo_expr = 1 # 
                          #auto_param_selection = T #
                          #verbose=F 
)

plot3 <- plot_cell_trajectory(mycds2, 
                              color_by = "orig.ident",
                              x=1,
                              show_tree = T,
                              show_backbone = F,
                              show_branch_points = T)
plot3

pdf(file="2F.pdf",width = 6,height = 4.5, onefile = FALSE)
print(plot3)
dev.off()