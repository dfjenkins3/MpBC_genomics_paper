rm(list=ls())
setwd("/restricted/projectnb/pathsig/work/yuqingz/emtSig_v2/")


########  Load Signature Data  ########
data_path <- "/restricted/projectnb/pathsig/work/dfj/20161207_emtSignatures/rsubread/"
dat1 <- read.table(paste(data_path, "emtSignatures_GFP_BCL2L11_SNAI_Twist.tpmlog", sep=""),
                   header=TRUE, row.names="Gene")
dat1 <- dat1[, -which(colnames(dat1)=="BCL2L11.9")]  # remove BCL2L11-9
dat2 <- read.table(paste(data_path, "emtSignatures_GFP_BECN.tpmlog", sep=""),
                   header=TRUE, row.names="Gene")
sum(rownames(dat1)!=rownames(dat2))

# combine two batches
sig_dat <- cbind(dat1, dat2)
batch_sig <- c(rep(1, ncol(dat1)), rep(2, ncol(dat2)))
# condition
geneLst <- c("GFP", "BCL2L11", "SNAI", "Twist", "BECN")
condition_sig <- rep(0, ncol(sig_dat))
for(i in 0:(length(geneLst)-1)){
  condition_sig[grep(geneLst[i+1], colnames(sig_dat))] <- i
}
table(condition_sig)

# indicator of genes that have 0 values in all samples within a batch
keep1 <- rowSums(abs(sig_dat[, batch_sig==1]))!=0
keep2 <- rowSums(abs(sig_dat[, batch_sig==2]))!=0
sig_dat <- sig_dat[keep1 & keep2, ]


# regular ComBat within signature
source("scripts/ComBat.R")
source("scripts/helper.R")
sig_dat_adj <- ComBat(dat=sig_dat, batch=batch_sig, 
                      mod=model.matrix(~as.factor(condition_sig)))


## load test set
load("new_val/val_Taube.RData")

overlap_genes <- intersect(rownames(sig_dat_adj), rownames(new.dat.mat))
sig_dat_adj <- sig_dat_adj[overlap_genes, ]
val_mat <- new.dat.mat[overlap_genes, ]
print(sum(rownames(sig_dat_adj)!=rownames(val_mat)))

## take over-express test sub-set
# change sample names to meaningful labels
print(colnames(val_mat) == meta.info$Sample)
colnames(val_mat) <- paste(meta.info$Pathway, meta.info$Rep, sep=".")
# take subset
sub_val_mat <- val_mat[, c(grep("pWZL", colnames(val_mat)),
                           grep("Snail",colnames(val_mat)),
                           grep("Twist",colnames(val_mat)))]


## ref ComBat between train & test, signature data as the reference
batch_idx <- c(rep(1, ncol(sig_dat_adj)), rep(2, ncol(sub_val_mat)))
adjdat_cmb <- ComBat(cbind(sig_dat_adj, sub_val_mat), batch=batch_idx, mod=NULL, ref.batch=1)
print(sum(adjdat_cmb[,1:ncol(sig_dat_adj)] != sig_dat_adj))


## Test sets 
# create input data list, each pathway as an element
sig_datLst <- list(BCL2L11=sig_dat_adj[, grep("BCL2L11", colnames(sig_dat_adj))],
                   BECN=sig_dat_adj[, grep("BECN", colnames(sig_dat_adj))],
                   SNAI=sig_dat_adj[, grep("SNAI", colnames(sig_dat_adj))],
                   Twist=sig_dat_adj[, grep("Twist", colnames(sig_dat_adj))],
                   GFP1=sig_dat_adj[, 1:10],
                   GFP2=sig_dat_adj[, 40:49])
val_indata <- append(sig_datLst, list(test=adjdat_cmb[, (ncol(sig_dat_adj)+1):ncol(adjdat_cmb)]))
print(sapply(val_indata, colnames))


########  ASSIGN  ########
library(ASSIGN)
ngene_seq <- c(10, 25, seq(from=0, to=500, by=50)[-1])
source("scripts/runassign.R")

excl_genes <- read.csv("freq_genes.csv", as.is=T)[, 1]
exclGeneLst <- list(BCL2L11=excl_genes, BECN=excl_genes, SNAI=excl_genes, Twist=excl_genes)

if(!dir.exists("./new_val/Excl_adST_2/")){
  dir.create("./new_val/Excl_adST_2/")
  dir.create("./new_val/Excl_adST_2/Val_Taube/")
}

for(ngene in ngene_seq){
  runassign(indata=val_indata, ngene=ngene, testname="Val_Taube",
            outputpathname="./new_val/Excl_adST_2/", exclude_geneList=exclGeneLst)
}  

