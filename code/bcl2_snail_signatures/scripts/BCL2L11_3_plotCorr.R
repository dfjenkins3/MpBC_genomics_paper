rm(list=ls())
setwd("/restricted/projectnb/pathsig/work/yuqingz/emtSig_v2/new_val/Incl_adST_2/Val_Elmore/")
load("../../val_Elmore.RData")

library(reshape)
library(ggplot2)

print(colnames(new.dat.mat) == meta.info$Sample)


####  BCL2 pathway 
BCL2_dirs <- dir()[grep("BCL2L11", dir())]
BCL2_dirs <- BCL2_dirs[order(sapply(strsplit(BCL2_dirs, "_"), function(strlst){return(as.numeric(strlst[2]))}))]
ngene_seq <- c(10, 25, seq(from=0, to=500, by=50)[-1])
# ngene_seq == sort(sapply(strsplit(BCL2_dirs, "_"), function(strlst){return(as.numeric(strlst[2]))}))


## read in predictions
BCL2_pred_res <- matrix(0, nrow=ncol(new.dat.mat), ncol=length(ngene_seq),
                        dimnames=list(colnames(new.dat.mat), ngene_seq))
for(i in seq_along(BCL2_dirs)){
  tmp <- read.csv(paste(BCL2_dirs[i], "pathway_activity_testset.csv", sep="/"), row.names="X")
  print(rownames(tmp) == rownames(BCL2_pred_res))
  BCL2_pred_res[, i] <- tmp[, "BCL2L11"]
  print(colnames(BCL2_pred_res)[i] == strsplit(BCL2_dirs[i], "_")[[1]][2]) # sanity check
}
write.csv(BCL2_pred_res, file="../Elmore_BCL2sig_predscores.csv", quote=FALSE)


## visualize
print(rownames(BCL2_pred_res) == meta.info$Sample)

# plot prediction
dat <- data.frame(Samples=rownames(BCL2_pred_res), BCL2_pred_res)
colnames(dat)[2:13] <- ngene_seq
dat2 <- melt(dat, id.vars="Samples")
colnames(dat2)[2:3] <- c("Signature_length", "Correlations")
dat3 <- merge(dat2, meta.info, by.x="Samples", by.y="Samples")

png("../Elmore_pred.png")
ggplot(dat3, aes(x=Group_abbv, y=Correlations)) + 
  geom_point() +
  facet_wrap(~Signature_length) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Predictions") + 
  ggtitle("BCL2 pathway predictions: Elmore data") +
  theme(axis.text.x = element_text(angle=60, hjust = 1))
dev.off()


# correlation with EC50
if(!dir.exists("../Elmore_corr_plots")){dir.create("../Elmore_corr_plots")}
corr_seq_pearson <- corr_seq_spearman <- c()
for(i in 1:ncol(BCL2_pred_res)){
  png(paste0("../Elmore_corr_plots/", colnames(BCL2_pred_res)[i], ".png"))
  plot(BCL2_pred_res[, i], meta.info$EC50.mean,
       xlab="Prediction", ylab="EC50", main=paste0(colnames(BCL2_pred_res)[i], "-gene Signature"))
  dev.off()
  corr_seq_pearson[i] <- cor(BCL2_pred_res[, i], meta.info$EC50.mean, method="pearson")
  corr_seq_spearman[i] <- cor(BCL2_pred_res[, i], meta.info$EC50.mean, method="spearman")
}
names(corr_seq_pearson) <- names(corr_seq_spearman) <- colnames(BCL2_pred_res)

dat <- data.frame(index=1:length(corr_seq_spearman), corr=corr_seq_spearman, 
                  names=names(corr_seq_spearman))
png("../Elmore_corr_EC50.png")
ggplot(data=dat, aes(x=index, y=corr)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(name="", labels=as.character(ngene_seq), breaks=1:12) +
  scale_y_continuous(name="Spearman Correlations") + 
  ggtitle("Correlation of BCL2 pathway prediction with EC50: Elmore data") +
  theme(axis.text=element_text(size=12))
dev.off()
  