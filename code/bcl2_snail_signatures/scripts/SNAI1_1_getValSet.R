rm(list=ls())
setwd("C:/Users/zhang/Google Drive/Johnson_Pathway/emt_signatures_reproduce/")
sapply(c("GEOquery", "hthgu133a.db"), require, character.only=TRUE)
# R -v 3.4.0, GEOquery 2.46.15, hthgu133a.db 3.2.3

#### Read data 
esets <- getGEO('GSE24202', GSEMatrix=TRUE)[[1]] 
dat.mat <- exprs(esets)
meta.info <- read.csv("GSE24202_Taube_meta.csv", header=TRUE)  


#### Log transform
dat.mat <- log2(dat.mat)
dat.mat[1:5, 1:5]


#### Annotation: match probe IDs with gene symbols
x <- hthgu133aSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])
head(xx)

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
identical(as.character(tmp[rownames(dat.mat[row_ind, ])]), symbol_vec)

# new data  
new.dat.mat <- dat.mat[row_ind, ]
rownames(new.dat.mat) <- symbol_vec


#### Save data
save(new.dat.mat, meta.info, file="val_Taube.RData")
