rm(list=ls())
setwd("/restricted/projectnb/pathsig/work/yuqingz/emtSig_v2/new_val/")


#### Read data 
library(GEOquery)
library(Biobase)
library(annotate)
library(hthgu133a.db)
esets <- getGEO('GSE24202', GSEMatrix=TRUE)[[1]] 
dat.mat <- exprs(esets)
meta.info <- read.csv("val_Taube_label.csv", header=TRUE)  

#library(oligo)
#library(pd.ht.hg.u133a)

#celfiles <- list.celfiles('./', listGzipped=TRUE)
#rawdata <- oligo::read.celfiles(celfiles)
#esets <- oligo::rma(rawdata)
#dat.mat <- exprs(esets)
#meta.info <- read.csv("val_Taube_label.csv", header=TRUE)  


#### Log transform??
dat.mat <- log2(dat.mat)
dat.mat[1:5, 1:5]


#### Annotation: match probe IDs with gene symbols
x <- hthgu133aSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])
head(xx)
#length(intersect(xx[,2], rownames(sig_dat)))

row_ind <- symbol_vec <- c()
for(i in 1:nrow(dat.mat)){
  ind <- match(rownames(dat.mat)[i], xx$probe_id)
  if(!is.na(ind)){
    row_ind <- c(row_ind, i)
    symbol_vec <- c(symbol_vec, xx$symbol[ind])
  }
}
# sanity check
tmp <- xx$symbol
names(tmp) <- xx$probe_id
sum(tmp[rownames(dat.mat[row_ind, ])] != symbol_vec)

## new data  
new.dat.mat <- dat.mat[row_ind, ]
rownames(new.dat.mat) <- symbol_vec


#### Save data
save(new.dat.mat, meta.info, file="val_Taube.RData")

