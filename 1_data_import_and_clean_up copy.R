# These RNA microarray analysis workflow contains these couple of sections: 
#"1_data_import_and_clean_up","2_QC","3_DGE","" 
## raw data are in the forms of CEL raw intensity files, this script imports all the files, do basic visualization, normalization, and import metadata sheet. 


# Data import -------------------------------------------------------------

library(oligo)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(tidyr)

# reading in cel files, rename them to shorter names
celFiles <- list.celfiles(path = "./PRO100776_UMich776_SEP/Raw_Data/", full.names = TRUE)
rawdata <- read.celfiles(celFiles)

# Some probe level visualizations: see how the data looks like 
head(exprs(rawdata), n = 2)
hist(rawdata)
boxplot(rawdata)

#sample metadata sheet
sdata <- read_xlsx("./PRO100776_UMich776_SEP/PRO100776_UMich776_SEP_Sample_Table")


# Normalization -----------------------------------------------------------
rawdata_norm <- rma(rawdata)
exprSet <- exprs(rawdata_norm)


# exprSet cleanup ---------------------------------------------------------

## annotate the probeset with gene symbols 

method <- "affy"
if (method == "affy") {
  #BiocManager::install("clariomshumanhttranscriptcluster.db")
  totalprbst <- dim(exprSet)[1]
  library(clariomshumanhttranscriptcluster.db)
  ids=toTable(clariomshumanhttranscriptclusterSYMBOL)
  exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,] #18599 83
  print(paste("annotated probeset no. / total probeset no. :" ,dim(exprSet)[1], "/", totalprbst))
  
  ids=ids[match(rownames(exprSet),ids$probe_id),] # match returns a vector of the positions of (first) matches of its first argument in its second.
  tmp = by(exprSet,ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] ) 
  # when multiple probe match to one gene, take the rowmeans and pick the maximum signal among all samples. 
  probes = as.character(tmp)
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  print(paste("probeset no. after clean-up :" ,dim(exprSet)[1]))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
}


## re-name the data 

sdata <- sdata[order(sdata$Well),]
names <- cbind(sdata$`Sample Name`)
colnames <- colnames(rawdata)
colnames <- sub("*.CEL", "", colnames)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
colnames <- substrRight(colnames, 3)
exprSet <- exprSet[,order(colnames)]
aligntest <- colnames[order(colnames)] == sdata$Well

if (FALSE %in% aligntest) {
  print("data not aligned, recheck ordering")
}

colnames(exprSet) <- names

## (optional) drop the irrelevant data (HeLa cell line was added by the processing facility as a loading control)

exprSet <- exprSet[ ,colnames(exprSet) != "HeLa"]
sdata <- sdata[sdata$`Sample Name` != "HeLa",]

## (optional) parse the condition column in sdata

sdata <- sdata %>% tidyr::separate(Condition, c("Type","Status"), ", ")

write.csv(exprSet,"./Organized Scripts/exprSet_step1.csv")
write.csv(sdata, "./Organized Scripts/sdata.csv")
