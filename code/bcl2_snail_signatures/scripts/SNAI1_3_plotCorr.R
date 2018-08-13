rm(list=ls())
setwd("/restricted/projectnb/pathsig/work/yuqingz/emtSig_v2/new_val/Excl_adST_2/Val_Taube/")
load("../../val_Taube.RData")
library(reshape)
library(ggplot2)

print(colnames(new.dat.mat) == meta.info$Sample)
colnames(new.dat.mat) <- paste(meta.info$Pathway, meta.info$Rep, sep=".")
# take subset
new.dat.mat <- new.dat.mat[, c(grep("pWZL", colnames(new.dat.mat)),
                               grep("Snail",colnames(new.dat.mat)),
                               grep("Twist",colnames(new.dat.mat)))]



####  SNAI pathway 
SNAI_dirs <- dir()[grep("SNAI", dir())]
SNAI_dirs <- SNAI_dirs[order(sapply(strsplit(SNAI_dirs, "_"), function(strlst){return(as.numeric(strlst[2]))}))]
ngene_seq <- c(10, 25, seq(from=0, to=500, by=50)[-1])
# ngene_seq == sort(sapply(strsplit(SNAI_dirs, "_"), function(strlst){return(as.numeric(strlst[2]))}))

## read in predictions
SNAI_pred_res <- matrix(0, nrow=ncol(new.dat.mat), ncol=length(ngene_seq),
                        dimnames=list(colnames(new.dat.mat), ngene_seq))
for(i in seq_along(SNAI_dirs)){
  tmp <- read.csv(paste(SNAI_dirs[i], "pathway_activity_testset.csv", sep="/"), row.names="X")
  print(rownames(tmp) == rownames(SNAI_pred_res))
  SNAI_pred_res[, i] <- tmp[, "SNAI"]
  print(colnames(SNAI_pred_res)[i] == strsplit(SNAI_dirs[i], "_")[[1]][2]) # sanity check
}

## visualize
#SNAI_pred_res <- cbind(meta.info, SNAI_pred_res)
#test <- melt(SNAI_pred_res)
# Bar plot
path_type <- c("pWZL", "Snail")
pltmat <- matrix(0, nrow=length(path_type), ncol=length(ngene_seq),
                 dimnames=list(path_type, ngene_seq))
for(i in seq_along(path_type)){
  tmp <- SNAI_pred_res[grep(path_type[i], rownames(SNAI_pred_res)), ]
  pltmat[path_type[i], ] <- colMeans(tmp)
}

pltmat <- melt(pltmat)
colnames(pltmat) <- c("Pathway", "N_gene", "Prediction")
pltmat[, "Pathway"] <- factor(as.character(pltmat[, "Pathway"]), levels=path_type)

png("../val_Taube_SNAI_2.png", width=700, height=700)
ggplot(pltmat, aes(y=Prediction, x=Pathway, color=Pathway, fill=Pathway), cex=2) + 
  geom_bar(stat="identity") +    
  facet_wrap(~N_gene) +
  #ylim(0, 0.4) +
  labs(title="Snail Pathway", x="", y="Pathway Activity Prediction", cex=2)
dev.off()

# Boxplot
tmp1 <- melt(SNAI_pred_res[1:6, ])
pltmat2 <- cbind(Pathway = factor(rep(path_type, each=3), levels=path_type),
                 Rep = as.factor(rep(c("rep1", "rep2", "rep3"), length(ngene_seq))),
                 tmp1)
colnames(pltmat2)[3:5] <- c("Sample", "N_gene", "Prediction")
print(pltmat2[,5]==tmp1[,3])
print(pltmat2[,4]==tmp1[,2])

png("../Box_Taube_SNAI_2.png", width=700, height=700)
p_snai <- ggplot(data=pltmat2, aes(x = Pathway, y = Prediction), fill=Pathway) +
  geom_boxplot() +
  facet_wrap(~N_gene) +
  scale_y_continuous(name = "Pathway Activity Prediction") + #,
                     #breaks = seq(0, 0.45, 0.05),
                     #limits=c(0, 0.45)) +
  scale_x_discrete(name = "") +
  ggtitle("Snail Pathway")
p_snai
dev.off()

write.csv(SNAI_pred_res, file="../Taube_SNAIsig_predscores.csv", quote=FALSE)


####  Twist pathway 
rm(list=ls())
load("../../val_Taube.RData")

print(colnames(new.dat.mat) == meta.info$Sample)
colnames(new.dat.mat) <- paste(meta.info$Pathway, meta.info$Rep, sep=".")
# take subset
new.dat.mat <- new.dat.mat[, c(grep("pWZL", colnames(new.dat.mat)),
                               grep("Snail",colnames(new.dat.mat)),
                               grep("Twist",colnames(new.dat.mat)))]


Twist_dirs <- dir()[grep("Twist", dir())]
Twist_dirs <- Twist_dirs[order(sapply(strsplit(Twist_dirs, "_"), function(strlst){return(as.numeric(strlst[2]))}))]
ngene_seq <- c(10, 25, seq(from=0, to=500, by=50)[-1])
# ngene_seq == sort(sapply(strsplit(Twist_dirs, "_"), function(strlst){return(as.numeric(strlst[2]))}))

## read in predictions
Twist_pred_res <- matrix(0, nrow=ncol(new.dat.mat), ncol=length(ngene_seq),
                         dimnames=list(colnames(new.dat.mat), ngene_seq))
for(i in seq_along(Twist_dirs)){
  tmp <- read.csv(paste(Twist_dirs[i], "pathway_activity_testset.csv", sep="/"), row.names="X")
  print(rownames(tmp) == rownames(Twist_pred_res))
  Twist_pred_res[, i] <- tmp[, "Twist"]
  print(colnames(Twist_pred_res)[i] == strsplit(Twist_dirs[i], "_")[[1]][2]) # sanity check
}

## visualize
#Twist_pred_res <- cbind(meta.info, Twist_pred_res)
#test <- melt(Twist_pred_res)
path_type <- c("pWZL", "Twist")
pltmat <- matrix(0, nrow=length(path_type), ncol=length(ngene_seq),
                 dimnames=list(path_type, ngene_seq))
for(i in seq_along(path_type)){
  tmp <- Twist_pred_res[grep(path_type[i], rownames(Twist_pred_res)), ]
  pltmat[path_type[i], ] <- colMeans(tmp)
}

pltmat <- melt(pltmat)
colnames(pltmat) <- c("Pathway", "N_gene", "Prediction")
pltmat[, "Pathway"] <- factor(as.character(pltmat[, "Pathway"]), levels=path_type)

png("../val_Taube_TWIST_2.png", width=700, height=700)
ggplot(pltmat, aes(y=Prediction, x=Pathway, color=Pathway, fill=Pathway), cex=2) + 
  geom_bar(stat="identity") +    
  facet_wrap(~N_gene) +
# ylim(0,0.35) +
  labs(title="TWIST Pathway", x="", y="Pathway Activity Prediction", cex=2)
dev.off()

# Boxplot
tmp2 <- melt(Twist_pred_res[c(1:3, 7:9), ])
pltmat2 <- cbind(Pathway = factor(rep(path_type, each=3), levels=path_type),
                 Rep = as.factor(rep(c("rep1", "rep2", "rep3"), length(ngene_seq))),
                 tmp2)
colnames(pltmat2)[3:5] <- c("Sample", "N_gene", "Prediction")
print(pltmat2[,5]==tmp2[,3])
print(pltmat2[,4]==tmp2[,2])

png("../Box_Taube_TWIST_2.png", width=700, height=700)
p_twist <- ggplot(data=pltmat2, aes(x = Pathway, y = Prediction), fill=Pathway) +
  geom_boxplot() +
  facet_wrap(~N_gene) +
  scale_y_continuous(name = "Pathway Activity Prediction") + #,
  #breaks = seq(0, 0.45, 0.05),
  #limits=c(0, 0.45)) +
  scale_x_discrete(name = "") +
  ggtitle("Twist Pathway")
p_twist
dev.off()

write.csv(Twist_pred_res, file="../Taube_TWISTsig_predscores.csv", quote=FALSE)
