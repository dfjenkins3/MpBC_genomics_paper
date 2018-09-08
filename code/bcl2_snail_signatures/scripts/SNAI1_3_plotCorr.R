rm(list=ls())
setwd("/restricted/projectnb/pathsig/work/yuqingz/emtSig_reproduce/results/Val_Taube/")
load("../../val_Taube.RData")
sapply(c("reshape", "ggplot2"), require, character.only=TRUE)  
#print(identical(colnames(new.dat.mat), as.character(meta.info$Sample)))

colnames(new.dat.mat) <- paste(meta.info$Pathway, meta.info$Rep, sep=".")
# take subset
new.dat.mat <- new.dat.mat[, c(grep("pWZL", colnames(new.dat.mat)),
                               grep("Snail",colnames(new.dat.mat)),
                               grep("Twist",colnames(new.dat.mat)))]

SNAI_dirs <- dir()[grep("SNAI", dir())]
SNAI_dirs <- SNAI_dirs[order(sapply(strsplit(SNAI_dirs, "_"), function(strlst){return(as.numeric(strlst[2]))}))]
ngene_seq <- c(10, 25, seq(from=0, to=500, by=50)[-1])


## read in predictions
SNAI_pred_res <- matrix(0, nrow=ncol(new.dat.mat), ncol=length(ngene_seq),
                        dimnames=list(colnames(new.dat.mat), ngene_seq))
for(i in seq_along(SNAI_dirs)){
  tmp <- read.csv(paste(SNAI_dirs[i], "pathway_activity_testset.csv", sep="/"), row.names="X")
  #print(identical(rownames(tmp), rownames(SNAI_pred_res)))
  SNAI_pred_res[, i] <- tmp[, "SNAI"]
  #print(identical(colnames(SNAI_pred_res)[i], strsplit(SNAI_dirs[i], "_")[[1]][2])) # sanity check
}


## visualize - Boxplot
path_type <- c("pWZL", "Snail")
pltmat2 <- cbind(Pathway = factor(rep(path_type, each=3), levels=path_type),
                 Rep = as.factor(rep(c("rep1", "rep2", "rep3"), length(ngene_seq))),
                 melt(SNAI_pred_res[1:6, ]))
colnames(pltmat2)[3:5] <- c("Sample", "N_gene", "Prediction")
#print(identical(pltmat2[,5], tmp1[,3])); print(identical(pltmat2[,4], tmp1[,2]))

png("../SNAI_validation_Taube.png", width=700, height=700)
ggplot(data=pltmat2, aes(x = Pathway, y = Prediction), fill=Pathway) +
  geom_boxplot() +
  facet_wrap(~N_gene) +
  scale_y_continuous(name = "Pathway Activity Prediction") + 
  scale_x_discrete(name = "") +
  ggtitle("Snail Pathway")
dev.off()

write.csv(SNAI_pred_res, file="../SNAIsig_predscores.csv", quote=FALSE)
