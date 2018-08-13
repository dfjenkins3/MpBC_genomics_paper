rm(list=ls())
setwd("/restricted/projectnb/pathsig/yuqingz/emtSig_v2/TCGA/")

#### Frequency of occurrence ####
dirLst <- dir()
geneLst <- c()

for(drname in dirLst){
  sigGenes <- read.csv(paste("./", drname, "/signature_gene_list_prior.csv", sep=""), as.is=T)
  geneLst <- c(geneLst, sigGenes[, 1])
}

# check geneLst
ngene_seq <- c(10, 25, seq(from=0, to=500, by=50)[-1])
sum(ngene_seq) * 4 == length(geneLst)

# find high freq genes
sorted.table <- sort(table(geneLst), decreasing=TRUE)
sublst <- sorted.table[sorted.table > (4*length(ngene_seq)*0.6)] 
freqGenes <- names(sublst)



#### Genes that appear in all pathways ####
sglpathway <- split(geneLst, ceiling(seq_along(geneLst)/sum(ngene_seq)))
# check sglpathway
sapply(sglpathway, length)
sum(sglpathway[[2]] != geneLst[(2785*1+1):(2785*2)])
TESTLst <- dirLst[grep("Twist", dirLst)]
testLst <- c()
for(drname in TESTLst){
  tmp <- read.csv(paste("./", drname, "/signature_gene_list_prior.csv", sep=""), as.is=T)
  testLst <- c(testLst, tmp[, 1])
}
sum(testLst != sglpathway[[4]])


sglpathTb <- lapply(sglpathway, unique)
# another check
length(unique(do.call(c, sglpathTb))) == length(unique(geneLst)) 

overlapGenes <- Reduce(intersect, sglpathTb)


#### Intersect between the two
noiseGenes <- intersect(freqGenes, overlapGenes)

# output frequency too
res <- sorted.table[match(noiseGenes, names(sorted.table))]
res <- round(sort(res, decreasing=T) / (4*length(ngene_seq)), 3)

# check
length(intersect(names(res), noiseGenes)) == length(noiseGenes)

output <- t(data.frame(as.list(res)))
write.csv(output, file="../freq_genes.csv")
