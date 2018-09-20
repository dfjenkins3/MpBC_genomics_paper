rm(list=ls())
setwd("/restricted/projectnb/pathsig/work/yuqingz/emtSig_reproduce/results/Val_Elmore/")
load("../../val_Elmore.RData")
sapply(c("reshape", "ggplot2"), require, character.only=TRUE)  
#print(identical(colnames(new.dat.mat), as.character(meta.info$Sample)))

BCL2_dirs <- dir()[grep("BCL2L11", dir())]
BCL2_dirs <- BCL2_dirs[order(sapply(strsplit(BCL2_dirs, "_"), function(strlst){return(as.numeric(strlst[2]))}))]
ngene_seq <- c(10, 25, seq(from=0, to=500, by=50)[-1])


## read in predictions
BCL2_pred_res <- matrix(0, nrow=ncol(new.dat.mat), ncol=length(ngene_seq),
                        dimnames=list(colnames(new.dat.mat), ngene_seq))
for(i in seq_along(BCL2_dirs)){
  tmp <- read.csv(paste(BCL2_dirs[i], "pathway_activity_testset.csv", sep="/"), row.names="X")
  #print(identical(rownames(tmp), rownames(BCL2_pred_res)))
  BCL2_pred_res[, i] <- tmp[, "BCL2L11"]
  #print(identical(colnames(BCL2_pred_res)[i], strsplit(BCL2_dirs[i], "_")[[1]][2])) # sanity check
}


## visualize
print(identical(rownames(BCL2_pred_res), as.character(meta.info$Sample)))
# correlation with EC50
corr_seq_pearson <- corr_seq_spearman <- c()
for(i in 1:ncol(BCL2_pred_res)){
  corr_seq_pearson[i] <- cor(BCL2_pred_res[, i], meta.info$EC50.mean, method="pearson")
  corr_seq_spearman[i] <- cor(BCL2_pred_res[, i], meta.info$EC50.mean, method="spearman")
}
names(corr_seq_pearson) <- names(corr_seq_spearman) <- colnames(BCL2_pred_res)

dat <- data.frame(index=1:length(corr_seq_spearman), corr=corr_seq_spearman, 
                  names=names(corr_seq_spearman))
png("../BCL2_validation_Elmore.png.png")
ggplot(data=dat, aes(x=index, y=corr)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(name="", labels=as.character(ngene_seq), breaks=1:12) +
  scale_y_continuous(name="Spearman Correlations") + 
  ggtitle("Correlation of BCL2 pathway prediction with EC50: Elmore data") +
  theme(axis.text=element_text(size=12))
dev.off()

write.csv(BCL2_pred_res, file="../BCL2sig_predscores.csv", quote=FALSE)
