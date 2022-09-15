#Part6
library(deconstructSigs)
library(maftools) 
library(ComplexHeatmap)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)
#################5A#################
load(file = 'maftools-LUAD.Rdata')  
#####TMB.txt#####
#准备突变频数
tmb = tmb(maf = laml)
tmb$total_perMB

TMB <- tmb[,c("Tumor_Sample_Barcode","total_perMB_log")]
colnames(TMB) <- c("ID","log10TMB")

write.table(TMB, "TMB.txt", row.names=F, col.names=T, sep="\t", quote=F)


#####mutsig.weightMatrix.txt"#####
maf <- as.data.frame(laml@data)
rmSilence = F 

if (rmSilence) {
  maf <- as.data.frame(maf[which(maf$Variant_Type == "SNP" & maf$Variant_Classification != "Silent"),]) #仅考虑SNP并移除沉默突变
} else {
  maf <- as.data.frame(maf[which(maf$Variant_Type == "SNP"),]) #
}
#maf$Chromosome <- paste0("chr",maf$Chromosome) # 求
maf$Tumor_Seq_Allele2
snp.count <- mut.to.sigs.input(mut.ref = maf, 
                               sample.id = "Tumor_Sample_Barcode", 
                               chr = "Chromosome", 
                               pos = "Start_Position", 
                               ref = "Reference_Allele", # 
                               alt = "Tumor_Seq_Allele2", # 
                               bsg = BSgenome.Hsapiens.UCSC.hg38) # 


cut.off <- 0.05 
mut.wt <- data.frame()
sigs.out.list <- list()
index <- 1
for (sample in rownames(snp.count)) {
  cat(paste0(sample," starts and ",length(rownames(snp.count))-index," samples remain to be analyzed!\n"))
  
  tmp <- whichSignatures(tumor.ref = snp.count, 
                         signatures.ref = signatures.cosmic, 
                         sample.id = sample, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome', 
                         signature.cutoff = cut.off)
  
  index <- index + 1
  #pdf(file.path(fig.path,paste0(sample,"_plotSignatures.pdf")))
  #plotSignatures(tmp)
  #invisible(dev.off())
  
  #pdf(file.path(fig.path,paste0(sample,"_weightPie.pdf")))
  #makePie(tmp)
  #invisible(dev.off())
  
  sigs.out.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights,unknown=tmp$unknown),row.names = sample)
  mut.wt <- rbind.data.frame(mut.wt,tmp)
}
write.table(mut.wt,"mutsig.weightMatrix.txt",sep = "\t",row.names = T,col.names = NA)


library(maftools)
library(stats)
library(tableone)

oncoplot(maf = laml,
         top = 1000,
         #genes = genelist,
         writeMatrix =T)

mut_matrix <- read.table("onco_matrix.txt",header = T,row.names = 1,sep = "\t",as.is = T)
mut_matrix[mut_matrix == "Missense_Mutation"] <- "1"
mut_matrix[mut_matrix == "Sense_Mutation"] <- "1"
mut_matrix[mut_matrix != "1"] <- "0"


mut <- mut_matrix
rownames(mut)
colnames(mut)
mut <- as.data.frame(t(mut))
mut$ID <- rownames(mut)

coln <- rownames(mut)
coln
coln <- gsub("[.]","-",coln)
rownames(mut) <- coln
rownames(mut) <- substr(rownames(mut),1,12)
mut$ID <- rownames(mut)

class <- subt
head(class)
class$ID <- rownames(class)
colnames(class) <- c("class","ID")
x <- merge(mut,class,by.x="ID",by.y="ID")

x <- as.data.frame(apply(x,2,as.character))
class(x$class)

var <- colnames(x[,2:1001])
result <- CreateTableOne(x, vars = var, strata = 'class')

result2 <- print(result)
#https://blog.csdn.net/baidu_33352045/article/details/80090214
write.csv(result2,"baseline.characteristics.CSV")

class(mut)
mydat <- as.data.frame(mut)
head(mydat)

mydat <- mydat[which(rownames(mydat)%in%class$ID),]
mydat <- merge(mydat,class,by.x="ID",by.y="ID")
rownames(mydat)
rownames(mydat) <- mydat$ID
mydat <- mydat[,-1]

coln <- rownames(mydat)
mydat <- as.data.frame(apply(mydat,2,as.numeric))
class(mydat$TP53)
mydat$class <- as.factor(mydat$class)
rownames(mydat) <- coln  

sum <- c()
mydata <-data.frame()

sum <- aggregate(mydat[,1], by=list(type = mydat$class),sum)
colnames(sum) <- c("type",colnames(mydat)[1])
mydata <- sum

for(i in 2:(ncol(mydat)-1)){
  sum <- aggregate(mydat[,i], by=list(type = mydat$class),sum)
  colnames(sum) <- c("type",colnames(mydat)[i])
  mydata <- merge(mydata,sum,by.x="type",by.y="type")
  print(i)
}

mydata[1:6,1:6]
xx <- as.data.frame(t(mydata))
colnames(xx) <- c("C1","C2","C3","C4")

head(xx)
xx <- xx[-1,]
head(xx)
coln <- rownames(xx)
xx <- as.data.frame(apply(xx,2,as.numeric))
rownames(xx) <- coln
xx$mut.percent <- (xx$C1 + xx$C2 + xx$C3 + xx$C4)/499
xx <- subset(xx, xx$mut.percent > 0.05)
selectgenes <- rownames(xx)

#xx <- xx[,-1]
#chisq <- chisq.test(xx[,2])
#p.result <- chisq$p.value
#X.squared <- chisq$statistic
#result.linshi<-cbind(i,X.squared,p.result)

result<-c()
xx <- mydata[which(colnames(mydata)%in%selectgenes)]

for (i in (1:ncol(xx))){
  newdat<-xx[,i]
  qq<-chisq.test(newdat)   #fisher.test()   chisq.test
  X.squared <- qq$statistic
  p.result <- qq$p.value
  result.linshi<-cbind(i,X.squared,p.result)
  result<-rbind(result,result.linshi)
}
gene<-names(xx)  
result<-cbind(gene,result)

result
write.csv(result,"mut.X.csv")


mut <- read.table("onco_matrix.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
colnames(mut) <- substr(colnames(mut),1,12)

mut[mut == 1] <- 1 #
mut[is.na(mut)] <- 0
mut[mut != 1] <- 0
#mut[mut == NA] <- 0 #
head(mut)

coln <- rownames(mut)
mut <- as.data.frame(apply(mut,2,as.character))
rownames(mut) <- coln
mut <- as.data.frame(t(mut))


tmb <- read.table("TMB.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(tmb) <- substr(rownames(tmb),1,12)
head(tmb)


mutsig <- read.table("mutsig.weightMatrix.txt",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(mutsig) <- substr(rownames(mutsig),1,12)
head(mutsig)


subt <- read.csv("scale50cell3.k=4.consensusClass.csv",
                 row.names = 1,header = F)
colnames(subt) <- c("CMOIC")
head(subt)


clust.col <- c("#E64B35","#00A087","#4DBBD5","#3C5488")
blue   <- "#5bc0eb"
red    <- "#f25f5c"


mutsig <- mutsig[,c("Signature.1","Signature.2","Signature.5","Signature.13")] # 
mutsig$APOBEC <- mutsig$Signature.2 + mutsig$Signature.13 # 
mutsig$CMOIC <- subt[rownames(mutsig),"CMOIC"] # 
mutsig <- mutsig[order(mutsig$CMOIC,-mutsig$APOBEC,decreasing = F),] # 
head(mutsig)

mutgene <- c("TP53",
             "KRAS",
             "KEAP1",
             "COL11A1",
             "CACNA1E",
             "HERC2",
             "CDH18",
             "PCDHB8",
             "DIDO1",
             "MYCBP2"
)


mut <- read.table("onco_matrix.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
colnames(mut) <- substr(colnames(mut),1,12)

onco.input <- mut[mutgene,rownames(mutsig),]
onco.input[onco.input == "Missense_Mutation"] <- "1"
onco.input[onco.input == "Sense_Mutation"] <- "1"
onco.input[onco.input != "Mutated"] <- "0"

head(onco.input)

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
  },
  Mutated = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#A60000", col = "#A60000")) 
  }
)
col = c("Mutated" ="#A60000") #

my_ann <- subt[rownames(mutsig),,drop = F]
coln <- rownames(my_ann)
my_ann <- as.data.frame(my_ann[-c(500:565),])
rownames(my_ann) <- coln[1:499]
colnames(my_ann) <- c("class")
class(my_ann$class)
my_ann$class <- as.factor(my_ann$class)

my_annotation = HeatmapAnnotation(df = my_ann, 
                                  col = list(CMOIC = c("1" = clust.col[1],
                                                       "2" = clust.col[2],
                                                       "3" = clust.col[3],
                                                       "4" = clust.col[4])))
）
top_anno <- anno_barplot(as.numeric(tmb[rownames(mutsig)[which(rownames(mutsig)%in%rownames(my_ann))],"log10TMB"]),
                         border = FALSE,
                         gp = gpar(fill = "#3379B4",border =NA,lty="blank"), 
                         height = unit(2.5, "cm"))


tmp <- mutsig[,c("Signature.2","Signature.13","Signature.1","Signature.5")]
tmp$Others <- 1 - rowSums(tmp) 
tmp <- tmp[which(rownames(tmp) %in% rownames(my_ann)),]

top_anno2 <- anno_barplot(as.matrix(tmp),
                          border = FALSE,
                          gp = gpar(fill = c(brewer.pal(6,"Paired")[c(2,1,6,5)],"grey90"), 
                                    border = NA, # 
                                    lty = "blank"),
                          height = unit(2, "cm")) # 

tmp <- as.data.frame(t(mut[mutgene,rownames(mutsig),]))
class(tmp$TP53)

coln <- rownames(tmp)
tmp <- as.data.frame(apply(tmp,2,as.numeric))
rownames(tmp) <- coln

mut.order <- names(sort(colSums(tmp),decreasing = T)) #
tmp$CMOIC <- subt[rownames(tmp),"CMOIC"]
pct <- NULL #
#i <- "TP53"
for (i in mut.order) {
  tmp1 <- tmp[,c(i,"CMOIC")]
  tmp1 <- as.data.frame.array(table(tmp1[,1],tmp1$CMOIC))[2,]/sum(tmp1[,1])
  pct <- rbind.data.frame(pct,tmp1)
}
rownames(pct) <- mut.order

right_anno <- anno_barplot(as.matrix(pct),
                           which = "row",
                           border = FALSE,
                           gp = gpar(fill = clust.col,border=NA,lty="blank"), 
                           bar_width = 0.6,
                           width = unit(1.8, "cm"),
                           height = unit(1, "cm"))

alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),
  Missense_Mutation = alter_graphic("rect", fill = "#E64B35FF"),
  Nonsense_Mutation = alter_graphic("rect", fill = mycolor[2]),
  Multi_Hit = alter_graphic("rect", fill = mycolor[3]),
  Splice_Site = alter_graphic("rect", fill = mycolor[4]),
  Frame_Shift_Del = alter_graphic("rect", fill = mycolor[5]),
  Frame_Shift_Ins = alter_graphic("rect", fill = mycolor[6]),
  In_Frame_Del = alter_graphic("rect", fill = mycolor[7]),
  In_Frame_Ins = alter_graphic("rect", height = 0.33, fill = mycolor[8])
)
library(ComplexHeatmap)
op1 <- oncoPrint(onco.input[mut.order,rownames(my_ann)], 
                 alter_fun = alter_fun,  
                 col = mycolor, #
                 bottom_annotation = NULL, 
                 top_annotation = c(HeatmapAnnotation(TMB = top_anno), 
                                    my_annotation, 
                                    HeatmapAnnotation(MutSig = top_anno2)), 
                 column_order = rownames(my_ann), 
                 #right_annotation = rowAnnotation(PCT = right_anno), 
                 show_pct = T, 
                 column_title = "", 
                 show_heatmap_legend = T, 
                 column_split = my_ann$class, 
                 column_title_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8),
                 alter_fun_is_vectorized = T)
op1



mut <- read.table("onco_matrix.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
colnames(mut) <- substr(colnames(mut),1,12)

onco.input <- mut[mutgene,rownames(mutsig),]
mut.order <- mutgene

alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),
  Missense_Mutation = alter_graphic("rect", fill = "#E64B35FF"),
  Nonsense_Mutation = alter_graphic("rect", fill = mycolor[2]),
  Multi_Hit = alter_graphic("rect", fill = mycolor[3]),
  Splice_Site = alter_graphic("rect", fill = mycolor[4]),
  Frame_Shift_Del = alter_graphic("rect", fill = mycolor[5]),
  Frame_Shift_Ins = alter_graphic("rect", fill = mycolor[6]),
  In_Frame_Del = alter_graphic("rect", fill = mycolor[7]),
  In_Frame_Ins = alter_graphic("rect", height = 0.33, fill = mycolor[8])
)
library(ComplexHeatmap)
op1 <- oncoPrint(onco.input[mut.order,rownames(my_ann)], # 
                 alter_fun = alter_fun,  # 
                 col = mycolor, # 
                 bottom_annotation = NULL, 
                 top_annotation = c(HeatmapAnnotation(TMB = top_anno), 
                                    my_annotation, 
                                    HeatmapAnnotation(MutSig = top_anno2)), 
                 column_order = rownames(my_ann), 
                 #right_annotation = rowAnnotation(PCT = right_anno), 
                 show_pct = T, 
                 column_title = "", 
                 show_heatmap_legend = T, 
                 column_split = my_ann$class, 
                 column_title_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8),
                 alter_fun_is_vectorized = T)
op1



library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(pheatmap)
library(stringr)
class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                  header = F, row.names = 1,check.names = F)
head(class)
class$ID  <- rownames(class)
colnames(class) <- c("class","ID")

all_cnv <- read.table("allclass_heatmap.txt",
                      row.names = 1,header = T,check.names = F,stringsAsFactors = F)
head(all_cnv)

coln <- colnames(all_cnv)[4:(ncol(all_cnv))]
coln
coln <- gsub("[.]","-",coln)
colnames(all_cnv)[4:(ncol(all_cnv))]<- coln

all_cnv <- all_cnv[,colnames(all_cnv)[substr(colnames(all_cnv)[4:(ncol(all_cnv))],1,12) %in% class$ID]] # 选取要展示的拷贝数变异

X =all_cnv
X[1:3,1:6]

X$seqnames  <- str_replace(X$seqnames, "chr","")
X$seqnames  <- str_replace(X$seqnames, "X","23")
class(X$seqnames)
X$seqnames <- as.numeric(X$seqnames)
X <- X[order(X[,1],X[,2]),]


X[1:3,1:6]

D <- X[,4:ncol(X)]
D[1:3,1:6]

D = D[,-which(colnames(D) %in% c("TCGA-CC-A8HV-01A-11D-A35Y-01"))]

D <- as.data.frame(t(D))

D = D[-which(rownames(D) %in% c("TCGA-38-4625-01A-01D-1204-01",
                                "TCGA-38-4626-01A-01D-1204-01",
                                "TCGA-38-4627-01A-01D-1204-01",
                                "TCGA-44-2655-01A-01D-0944-01",
                                "TCGA-44-2656-01A-02D-0944-01",
                                "TCGA-44-2656-01A-02D-A273-01",
                                "TCGA-44-2659-01A-01D-0944-01",
                                "TCGA-44-2661-01A-01D-1549-01",
                                "TCGA-44-2662-01A-01D-0944-01",
                                "TCGA-44-2662-01A-01D-A273-01",
                                "TCGA-44-2665-01A-01D-A273-01",
                                "TCGA-44-2665-01A-01D-0944-01",
                                "TCGA-44-2668-01A-01D-0944-01",
                                "TCGA-44-2668-01A-01D-A273-01",
                                "TCGA-44-3396-01A-01D-1204-01",
                                "TCGA-44-3398-01A-01D-1549-01",
                                "TCGA-44-3918-01A-01D-A273-01",
                                "TCGA-44-4112-01A-01D-A273-01",
                                "TCGA-44-5645-01A-01D-A273-01",
                                "TCGA-44-6146-01A-11D-A273-01",
                                "TCGA-44-6147-01A-11D-A273-01",
                                "TCGA-44-6775-01A-11D-A273-01",
                                "TCGA-44-2666-01A-01D-A273-01",
                                "TCGA-44-2657-01A-01D-1549-01")),]

D[1:3,1:6]
rownames(D) <- substr(rownames(D),1,12)
D$ID <- rownames(D)
D <- merge(D,class,by.x = "ID",by.y = "ID")
class(D$class)

D$class <- as.numeric(D$class)
D <- D[order(D$class),]
rownames(D) <- D$ID
D = D[,-which(colnames(D) %in% c("ID","class"))]
#coln <- rownames(D)
D <- as.data.frame(t(D))


dim(D)
colnames(D)
rownames(D)

library(RColorBrewer)
#cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256) 

library(stringr)
ann_col <- class[(rownames(class)%in%colnames(D)),] #
annotation_col = as.data.frame(ann_col[,1])
rownames(annotation_col) =  rownames(ann_col)
colnames(annotation_col) <- c("class")

annotation_row = data.frame(
  chr= factor(X$seqnames, levels = unique(X$seqnames))
)
rownames(annotation_row) = rownames(X)

#annotation_row <- annotation_row[order(annotation_row[,1]),]
dim(D)
dim(annotation_col)
dim(annotation_row)

D <- as.data.frame(D)
D <- as.matrix(D)

range(D)
D[1:6,1:6]

library(circlize)
col_fun <- colorRamp2(
  c(-1.5, 0, 4), 
  c("#3C5488", "#F8F8FF", "#E64B35")
)

bk <- c(seq(-1.5,-0.1,by=0.01),seq(0,4,by=0.01))

table(X[,2])
setwd("F:\\frontiers\\TMB&CNV\\CNV")
p <- pheatmap(D,
              cluster_cols = F,cluster_rows = F,
              #col=rev(cols),
              annotation_col = annotation_col,
              annotation_row = annotation_row,
              show_rownames=F,
              show_colnames=F,
              gaps_col = c(78,223,390),
              gaps_row = c(2339,3818,5148,6026,7040,8193,9281,10118,11055,11970,13428,14592,15143,15936,16756,17719,19087,19440,21024,21697,22004,22536),
              #color = col_fun,
              color = c(colorRampPalette(colors = c("#3C5488","#F8F8FF"))(length(bk)/2),colorRampPalette(colors = c("#F8F8FF","#E64B35"))(length(bk)/2)),
              #breaks=bk,
              border=FALSE,
              use_raster = T,
              annotation_legend = T  #去注释图例
              
) 
p

pdf("5A.pdf",21,29.7)
print(p)
dev.off()

#####5B#####
class <- read.csv("scale50cell3.k=4.consensusClass.csv",
                  header = F, row.names = 1,check.names = F)
head(class)
class$ID  <- rownames(class)
colnames(class) <- c("class","ID")

tumorCNV <- read.table("segment_file.txt",sep = "\t",header = T)
head(tumorCNV)

class1_id <- class[class$class=="1", "ID"]
class2_id <- class[class$class=="2", "ID"]
class3_id <- class[class$class=="3", "ID"]
class4_id <- class[class$class=="4", "ID"]

# Get subtype segments information
tumorCNV$x <- substr(tumorCNV$Sample,1,12)
class1_seg <- tumorCNV[tumorCNV$x %in% class1_id,]
class2_seg <- tumorCNV[tumorCNV$x %in% class2_id,]
class3_seg <- tumorCNV[tumorCNV$x %in% class3_id,]
class4_seg <- tumorCNV[tumorCNV$x %in% class4_id,]

head(class1_seg)

class1_seg <- class1_seg[,-7]
class2_seg <- class2_seg[,-7]
class3_seg <- class3_seg[,-7]
class4_seg <- class4_seg[,-7]

head(class1_seg)

write.table(class1_seg, file="class1_seg.txt", sep="\t", row.names=F, quote = F)
write.table(class2_seg, file="class2_seg.txt", sep="\t", row.names=F, quote = F)
write.table(class3_seg, file="class3_seg.txt", sep="\t", row.names=F, quote = F)
write.table(class4_seg, file="class4_seg.txt", sep="\t", row.names=F, quote = F)


#分组绘图
library(maftools)

lihc.gistic <- readGistic(gisticAllLesionsFile="all_lesions.conf_90.txt", 
                          gisticAmpGenesFile="amp_genes.conf_90.txt", 
                          gisticDelGenesFile="del_genes.conf_90.txt", 
                          gisticScoresFile="scores.gistic", isTCGA=TRUE)
getSampleSummary(lihc.gistic)
getGeneSummary(lihc.gistic)
getCytobandSummary(lihc.gistic)

write.GisticSummary(gistic=lihc.gistic, basename="luad_gistic2_class4")
pdf("class4CNV.fdrCutOff0.01.pdf",16,4.5)
gisticChromPlot(gistic=lihc.gistic, markBands="all",ref.build = "hg38",
                y_lims = c(1,-0.5),fdrCutOff = 0.01,txtSize=2
                #,	maf = maf,mutGenes = "TP53"
)
dev.off()






