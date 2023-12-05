# note that some of these QC uses the rawdata (type: ExpressionFeatureSet) in step1. 
# This script does the quality control of each samples. That include: NUSE, RLE, PCA, filtering

newDirectory <- paste(getwd(), "/","RDART_0331data_QC", sep="")
dir.create(newDirectory)

# NUSE,RLE ----------------------------------------------------------------

plmFit <- fitProbeLevelModel(rawdata)
n_cel_files <- length(celFiles)
png(file.path(newDirectory, "plmFit-NUSE-scaled.png"),width = 50 * n_cel_files, height = 800)
boxplot(plmFit, main = "NUSE", ylim = c(0.95,1.5))
dev.off()

png(file.path(newDirectory, "plmFit-RLE.png"),
    width = 50 * n_cel_files, height = 800)
RLE(plmFit, main = "RLE", cex.axis = 1.2, cex.lab = 1.2)
dev.off()


# PCA ---------------------------------------------------------------------
library(ggplot2)
library(gridExtra)

pdf(file.path(newDirectory,"PCA_all.pdf"))

## plot 1: cell vs exosome samples

sample_type <- TRUE
if(sample_type){
  tMat <- t(as.matrix(exprSet))
  Sacp <- prcomp(tMat)
  aa <- as.data.frame(Sacp$x)
  aa$group <- sdata$`Sample Type`
  p <- ggplot(aa,aes(x = PC1, y = PC2, color = group))
  p1 <- p + geom_point() + theme(legend.position = "top")
}

## plot 2: healthy, cell lines, patient. 

CTC_sample_type <- TRUE
if(CTC_sample_type){
  tMat <- t(as.matrix(exprSet))
  Sacp <- prcomp(tMat)
  aa <- as.data.frame(Sacp$x)
  aa$group <- sdata$Type
  p <- ggplot(aa,aes(x = PC1, y = PC2, color = group))
  p2 <- p + geom_point() + theme(legend.position = "top")
}

## plot 3: within only CTC samples and healthy, patient time point. (note 20210629 need to change this)

sdata_CTC <- sdata[sdata$`Sample Type` == "Cell",]
sdata_CTC <- sdata_CTC[sdata_CTC$Type %in% c("Tumor", "Healthy"),]

exprSet_CTC <- exprSet[ ,colnames(exprSet) %in% sdata_CTC$`Sample Name`]

timepoint <- TRUE
if(timepoint){
  tMat <- t(as.matrix(exprSet_CTC))
  Sacp <- prcomp(tMat)
  aa <- as.data.frame(Sacp$x)
  aa$group <- sdata_CTC$TimePoint
  p <- ggplot(aa,aes(x = PC1, y = PC2, color = group))
  p3 <- p + geom_point() + theme(legend.position = "top")
}

## plot 4: within only CTC samples and healthy, by patient status

sdata_CTC <- sdata[sdata$`Sample Type` == "Cell",]
sdata_CTC <- sdata_CTC[sdata_CTC$Type %in% c("Tumor", "Healthy"),]

exprSet_CTC <- exprSet[ ,colnames(exprSet) %in% sdata_CTC$`Sample Name`]

status <- TRUE
if(status){
  tMat <- t(as.matrix(exprSet_CTC))
  Sacp <- prcomp(tMat)
  aa <- as.data.frame(Sacp$x)
  aa$group <- sdata_CTC$Status
  p <- ggplot(aa,aes(x = PC1, y = PC2, color = group))
  p4 <- p + geom_point() + theme(legend.position = "top")
}

grid.arrange(p1, p2, p3, p4, ncol = 2)

graphics.off()



