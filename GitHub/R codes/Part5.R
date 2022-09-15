#Part5
#####差异分析，LIMMA#####

library(limma)

createList <- function(group=NULL) {
  
  tumorsam <- names(group)
  sampleList = list()
  treatsamList =list()
  treatnameList <- c()
  ctrlnameList <- c()
  
  #A-1: 
  sampleList[[1]] = tumorsam
  treatsamList[[1]] = intersect(tumorsam, names(group[group=="1"])) # 
  treatnameList[1] <- "1" # 
  ctrlnameList[1] <- "Others" # 
  
  #A-2:
  sampleList[[2]] = tumorsam
  treatsamList[[2]] = intersect(tumorsam, names(group[group=="2"]))
  treatnameList[2] <- "2"
  ctrlnameList[2] <- "Others"
  
  #A-3
  sampleList[[3]] = tumorsam
  treatsamList[[3]] = intersect(tumorsam, names(group[group=="3"]))
  treatnameList[3] <- "3"
  ctrlnameList[3] <- "Others"
  
  sampleList[[4]] = tumorsam
  treatsamList[[4]] = intersect(tumorsam, names(group[group=="4"]))
  treatnameList[4] <- "4"
  ctrlnameList[4] <- "Others"
  
  
  
  return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
}

DEGs_limma <- function(datas, groups, case, ref) {
  up_index <- which(groups[, 2] == case)
  groups[, 2][up_index] <- 'T0'
  down_index <- which(groups[, 2] == ref)
  groups[, 2][down_index] <- 'T1'
  
  groups <- groups[c(up_index, down_index), ] # 
  datas <- datas[, groups[,1]]           # 
  group_list <- factor(groups[, 2])
  design <- model.matrix(~group_list)
  fit <- lmFit(datas, design)
  fit <- eBayes(fit)
  options(digits = 4)
  topTable(fit, coef=2, adjust='BH')
  res <- topTable(fit, coef=2, adjust='BH', number = Inf)
  res$logFC <- -res$logFC
  return(res)
}

twoclasslimma <- function(res.path=NULL, countsTable=NULL, prefix=NULL, complist=NULL, overwt=FALSE) {
  
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(countsTable)
  
  options(warn=1)
  for (k in 1:length(sampleList)) { # 
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]] 
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    
    compname <- paste(treatname, "_vs_", ctrlname, sep="") # 
    tmp = rep("others", times=length(allsamples))
    names(tmp) <- allsamples
    tmp[samples]="control"
    tmp[treatsam]="treatment"
    outfile <- file.path( res.path, paste(prefix, "_limma_test_result.", compname, ".txt", sep="") )
    if (file.exists(outfile) & (overwt==FALSE)) { # 
      cat(k, ":", compname, "exists and skipped;\n")
      next
    }
    
    saminfo <- data.frame("Type"= tmp[samples],"SampleID"= samples,stringsAsFactors = F)
    cts <- countsTable[,samples]
    coldata <- saminfo[samples,]
    colnames(coldata) <- c("Groups","Samples")
    
    a <- rownames(cts)
    b <- colnames(cts)
    cts = as.data.frame(lapply(cts,as.numeric))
    rownames(cts) <- a
    colnames(cts) <- b
    # 
    dat_limma <- DEGs_limma(datas = cts,
                            groups = coldata[, c("Samples", "Groups")],
                            case = 'treatment', # 
                            ref = 'control')
    dat_limma$ID <- rownames(dat_limma)
    dat_limma <- dat_limma[,c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
    colnames(dat_limma) <- c("ID","log2FC","AveExpr","t","P.Value","FDR","B")
    
    #
    write.table(dat_limma, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
    cat(k, ",")
  }
  options(warn=0)
}


expr <- read.csv("F:\\肺癌第一次返修\\genematrixTCGA + GSE68465\\combined.expr.combat.csv",
                 sep = ",",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
range(expr)

subt <-  read.csv("scale50cell3.k=4.consensusClass.csv",
                  header = F, row.names = 1,check.names = F)
head(subt)
colnames(subt) <-  ("TCGA_Subtype")
n.sub.label <- unique(subt$TCGA_Subtype) 
n.sub.label
n.sub <- length(table(subt$TCGA_Subtype)) 
n.sub
expr <- expr[,which(colnames(expr) %in% rownames(subt))]


group <- subt$TCGA_Subtype
names(group) <- rownames(subt) 
complist <- createList(group=group)


twoclasslimma(res.path = ".", 
              countsTable = expr, #expr[index,intersect(colnames(expr),rownames(subt))],
              prefix = "LUAD", 
              complist = complist,
              overwt = F)


#####4B&C#####
library(GSVA)
library(clusterProfiler) 
library(pheatmap)
library(gplots)
DEfiles <- c("LUAD_limma_test_result.1_vs_Others.txt",
             "LUAD_limma_test_result.2_vs_Others.txt",
             "LUAD_limma_test_result.3_vs_Others.txt",
             "LUAD_limma_test_result.4_vs_Others.txt"
)
# 设置阈值
fdrcut <- 0.05 
#logfccut <- 2 # 理论上不需要这个阈值
top <- 1000 

degs.list <- list()
for (i in 1:n.sub) {
  degs <- read.table(DEfiles[i],sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
  head(degs)
  degs.list[[n.sub.label[i]]] <- as.data.frame(na.omit(degs))
}


subtype_specific_gsea <- function(msigdb=NULL,n.top=10,mode=c("up","down"),degs.list=NULL,subtype.label=NULL,nPerm.gsea=1000,minGSSize.gsea=10,maxGSSize.gsea=500,pvalueCutoff.gsea=1){
  
  MSigDB <- read.gmt(msigdb)
  GSEA.list <- top.gs <- list() #
  
  if(!is.element(mode, c("up", "dn"))) { stop("mode must be up or dn!\n") }
  
  for (i in 1:n.sub) {
    degs <- degs.list[[n.sub.label[i]]]
    geneList <- degs$log2FC; names(geneList) <- as.character(rownames(degs))
    geneList <- sort(geneList,decreasing = T) # ranked gene set
    # 
    cat(paste0("GSEA for ",n.sub.label[i]," starts!\n"))
    #GSEA.list[[subtype.label[i]]] 
    GSEA.list[[i]]<- GSEA(geneList = geneList,
                          TERM2GENE=MSigDB,
                          nPerm = nPerm.gsea,
                          minGSSize = minGSSize.gsea,
                          maxGSSize = maxGSSize.gsea,
                          seed = T,
                          verbose = F,
                          pvalueCutoff = pvalueCutoff.gsea) # 
    
    #GSEA.dat <- as.data.frame(GSEA.list[[subtype.label[i]]])
    GSEA.dat <- as.data.frame(GSEA.list[i])
    if(mode == "up") {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = T),] # 
    } else {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = F),] # 
    }
    
    # 
    write.table(GSEA.dat,paste0(n.sub.label[i],"_degs_",mode,"_gsea.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
    
    # 
    top.gs[[i]] <- rownames(GSEA.dat)[1:n.top] 
  }
  
  # 
  gs <- list()
  for (i in as.character(unlist(top.gs))) {
    gs[[i]] <- MSigDB[which(MSigDB$term %in% i),"gene"]
  }
  
  return(list(mode=mode,top.gs=top.gs,gs=gs))
}

jco <- c("#F2CCCC","#E6D8CF","#D5E3F0","#FDE7DA","#E2D6EC", "#CCEFDB")
annColors <- list(subtype=c("C1"=jco[1],"C2"=jco[2],"C3"=jco[3],"C4"=jco[4]))

msigdfFile = "h.all.v7.4.symbols.gmt"
n.top = 10

mode = "up" 

gs.up <- subtype_specific_gsea(msigdb = msigdfFile,
                               n.top = n.top,
                               degs.list = degs.list,
                               subtype.label = n.sub.label,
                               mode = mode)

gsva_gs.up <- gsva(as.matrix(expr), gs.up$gs, method="gsva") 
dim(gsva_gs.up)


gsva_gs.up_mean <- data.frame(row.names = rownames(gsva_gs.up)) 
for (i in n.sub.label) {
  gsva_gs.up_mean <- cbind.data.frame(gsva_gs.up_mean,
                                      data.frame(rowMeans(gsva_gs.up[,rownames(subt)[which(subt$TCGA_Subtype == i)]])))
}
colnames(gsva_gs.up_mean) <- n.sub.label

annRows <- data.frame(subtype=rep(n.sub.label,each=n.top), names = unlist(gs.up$top.gs), stringsAsFactors = F)
annRows <- annRows[!duplicated(annRows$names),]; rownames(annRows) <- annRows$names # 
filename <- paste0("HALL","_specific_top_",mode,"_gsea.pdf")
pheatmap(gsva_gs.up_mean[rownames(annRows),],
         cellwidth = 10, cellheight = 10,
         #color = bluered(64), #
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA, #
         annotation_row = annRows[,"subtype",drop = F],
         annotation_colors = annColors,#,
         filename = filename
)


mode = "dn"
gs.dn <- subtype_specific_gsea(msigdb = msigdfFile,
                               n.top = n.top,
                               degs.list = degs.list,
                               subtype.label = n.sub.label,
                               mode = mode)
# 
gsva_gs.dn <- gsva(as.matrix(expr), gs.dn$gs, method="gsva") # 
dim(gsva_gs.dn)

# 
gsva_gs.dn_mean <- data.frame(row.names = rownames(gsva_gs.dn)) 
for (i in n.sub.label) {
  gsva_gs.dn_mean <- cbind.data.frame(gsva_gs.dn_mean,
                                      data.frame(rowMeans(gsva_gs.dn[,rownames(subt)[which(subt$TCGA_Subtype == i)]])))
}
colnames(gsva_gs.dn_mean) <- n.sub.label

#
annRows <- data.frame(subtype=rep(n.sub.label,each=n.top), names = unlist(gs.dn$top.gs), stringsAsFactors = F)
annRows <- annRows[!duplicated(annRows$names),]; rownames(annRows) <- annRows$names # 倘若出现一条通路在>=2个亚型中下调，去掉重复值，这种情况在亚型较多的时候会发生

filename <- paste0("HALL","_specific_top_",mode,"_gsea.pdf")
pheatmap(gsva_gs.dn_mean[rownames(annRows),],
         cellwidth = 10, cellheight = 10,
         #color = bluered(64), #
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA, #
         annotation_row = annRows[,"subtype",drop = F],
         annotation_colors = annColors,
         filename = filename
)


#####4A#####
library(ggplot2)

#1
deg.data <- read.table('LUAD_limma_test_result.1_vs_Others.txt',sep="\t", header = T,stringsAsFactors = F, check.names = F)
outTab2 <- deg.data

outTab2$reg <- ifelse(outTab2$log2FC >1 &  outTab2$FDR <0.05, 'Up', ifelse(outTab2$log2FC < -1 &  outTab2$FDR <0.05, 'Down', 'NON'))

mycol <- c("#185AA8", "grey", "#B83D3C", "#8C510A") #cluster的颜色，如果有更多类，就给更多的颜色

p <- ggplot(data = outTab2, aes(x= log2FC, y = -log10(FDR) )) +
  geom_point(data = outTab2,aes(x= log2FC, y = -log10(FDR), color= factor(reg)),size = 3)+
  scale_color_manual(values = mycol)+ #ylim(0,3)+
  theme_bw()+ geom_vline( xintercept = c(-1,1)) + geom_hline(yintercept = -log10(0.05))+
  ggrepel::geom_text_repel(data = outTab2[which(outTab2$ID %in% c('SFTPC','CDC20',
                                                                  'C7','MYBL2',
                                                                  'ACKR1','NEK2',
                                                                  'SCGB1A1','S100P',
                                                                  'ADH1B','FOXM1'
  )),] ,size= 4,aes(label= ID ),
  box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"), 
  segment.color = "black", show.legend = FALSE)

pdf(file="1.pdf",width = 8,height = 4, onefile = FALSE)
print(p)
dev.off()


#2
deg.data <- read.table('LUAD_limma_test_result.2_vs_Others.txt',sep="\t", header = T,stringsAsFactors = F, check.names = F)
outTab2 <- deg.data

outTab2$reg <- ifelse(outTab2$log2FC >1 &  outTab2$FDR <0.05, 'Up', ifelse(outTab2$log2FC < -1 &  outTab2$FDR <0.05, 'Down', 'NON'))

p<- ggplot(data = outTab2, aes(x= log2FC, y = -log10(FDR) )) +
  geom_point(data = outTab2,aes(x= log2FC, y = -log10(FDR), color= factor(reg)),size = 3)+
  scale_color_manual(values = mycol)+ #ylim(0,3)+
  theme_bw()+ geom_vline( xintercept = c(-1,1)) + geom_hline(yintercept = -log10(0.05))+
  ggrepel::geom_text_repel(data = outTab2[which(outTab2$ID %in% c('SFTPC','FGB',
                                                                  'C7','FGA',
                                                                  'MFAP4','MAGEA6',
                                                                  'SCGB1A1','FGL1',
                                                                  'CHIT1','MYBL2'
  )),] ,size= 4,aes(label= ID ),
  box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"), 
  segment.color = "black", show.legend = FALSE)

pdf(file="2.pdf",width = 8,height = 4, onefile = FALSE)
print(p)
dev.off()

#3
deg.data <- read.table('LUAD_limma_test_result.3_vs_Others.txt',sep="\t", header = T,stringsAsFactors = F, check.names = F)
outTab2 <- deg.data

outTab2$reg <- ifelse(outTab2$log2FC >1 &  outTab2$FDR <0.05, 'Up', ifelse(outTab2$log2FC < -1 &  outTab2$FDR <0.05, 'Down', 'NON'))

p <- ggplot(data = outTab2, aes(x= log2FC, y = -log10(FDR) )) +
  geom_point(data = outTab2,aes(x= log2FC, y = -log10(FDR), color= factor(reg)),size = 3)+
  scale_color_manual(values = mycol)+ #ylim(0,3)+
  theme_bw()+ geom_vline( xintercept = c(-1,1)) + geom_hline(yintercept = -log10(0.05))+
  ggrepel::geom_text_repel(data = outTab2[which(outTab2$ID %in% c('CXCL13','C1orf116',
                                                                  'CXCL9','SFTPC',
                                                                  'NKG7','CYP4B1',
                                                                  'ADAMDEC1','CYP2B7P',
                                                                  'CXCL10','PGC'
  )),] ,size= 4,aes(label= ID ),
  box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"), 
  segment.color = "black", show.legend = FALSE)

pdf(file="3.pdf",width = 8,height = 4, onefile = FALSE)
print(p)
dev.off()

#4
deg.data <- read.table('LUAD_limma_test_result.4_vs_Others.txt',sep="\t", header = T,stringsAsFactors = F, check.names = F)
outTab2 <- deg.data

outTab2$reg <- ifelse(outTab2$log2FC >1 &  outTab2$FDR <0.05, 'Up', ifelse(outTab2$log2FC < -1 &  outTab2$FDR <0.05, 'Down', 'NON'))

p<- ggplot(data = outTab2, aes(x= log2FC, y = -log10(FDR) )) +
  geom_point(data = outTab2,aes(x= log2FC, y = -log10(FDR), color= factor(reg)),size = 3)+
  scale_color_manual(values = mycol)+ #ylim(0,3)+
  theme_bw()+ geom_vline( xintercept = c(-1,1)) + geom_hline(yintercept = -log10(0.05))+
  ggrepel::geom_text_repel(data = outTab2[which(outTab2$ID %in% c('CXCL13','CYP2B7P',
                                                                  'CXCL9','SFTPC',
                                                                  'NKG7','CYP4B1',
                                                                  'ADAMDEC1','CACNA2D2',
                                                                  'CXCL10','PGC'
  )),] ,size= 4,aes(label= ID ),
  box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"), 
  segment.color = "black", show.legend = FALSE)

pdf(file="4.pdf",width = 8,height = 4, onefile = FALSE)
print(p)
dev.off()
