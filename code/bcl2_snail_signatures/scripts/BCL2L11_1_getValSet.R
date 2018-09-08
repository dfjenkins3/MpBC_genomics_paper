rm(list=ls())
setwd("C:/Users/zhang/Google Drive/Johnson_Pathway/emt_signatures_reproduce/")
sapply(c("GEOquery", "hgu133a2.db"), require, character.only=TRUE)
# R -v 3.

## Load data 
# download expression set from GEO
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10841 # data location
esets <- getGEO('GSE10841')[[1]] 
dat.mat <- exprs(esets)

# read in sample information
meta.info <- read.csv("GSE10841_Elmore_meta.csv")
identical(colnames(dat.mat), as.character(meta.info$Sample))  # sample name match between meta and expression data

# remove NCI-H740 (EC50 >10, no clear result)
dat.mat <- dat.mat[, -grep(">10", meta.info[,4])]
meta.info <- meta.info[-grep(">10", meta.info[,4]), ]
meta.info[,4] <- as.numeric(as.character(meta.info[,4]))

# log transform
dat.mat <- log2(dat.mat)


## Annotation
x <- hgu133a2SYMBOL 
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])

row_ind <- symbol_vec <- c()
for(i in 1:nrow(dat.mat)){
  ind <- match(rownames(dat.mat)[i], xx$probe_id)
  if(!is.na(ind)){
    row_ind <- c(row_ind, i)
    symbol_vec <- c(symbol_vec, xx$symbol[ind])
  }
}

# sanity check of annotation
tmp <- xx$symbol
names(tmp) <- xx$probe_id
identical(as.character(tmp[rownames(dat.mat[row_ind, ])]), symbol_vec)


## Save annotated data
new.dat.mat <- dat.mat[row_ind, ]
rownames(new.dat.mat) <- symbol_vec

save(new.dat.mat, meta.info, file="val_Elmore.RData")
