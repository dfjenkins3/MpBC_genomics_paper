###install NanoStringDiff
source("https://bioconductor.org/biocLite.R")
biocLite("NanoStringDiff")

###specify directory and path to counts file
directory <- "C:/Users/jasmi/Documents/Bild Lab/NanoString/Yuan MpBC/R Code/"
path <- paste(directory,"metaplastic_tnbc_coh_labels.csv", sep="/")

###load the phenotype data
designs=read.csv(file="C:/Users/jasmi/Documents/Bild Lab/NanoString/Yuan MpBC/R Code/metaplastic_tnbc_phenotype_labels.csv", header=TRUE, row.names=1, sep=",")
designs

###Call NanoStringDiff and create the NS set
library("NanoStringDiff")
NanoStringData=createNanoStringSetFromCsv(path,header=TRUE,designs)

###Make a matrix to prepare for normalization, that only considers one phenotype label, in this case group
pheno=pData(NanoStringData)
group=pheno$group
design.full=model.matrix(~0+group)
design.full

###Normalization of the data
NanoStringData=estNormalizationFactors(NanoStringData)
positiveFactor(NanoStringData)
negativeFactor(NanoStringData)
housekeepingFactor(NanoStringData)

###perform glm ratio test of MpBC to control TNBC
###control group gets contrast -1 designation, group of interest gets +1
result=glm.LRT(NanoStringData,design.full,contrast=c(1,-1))
head(result$table)
str(result)

###write result out to csv
write.csv(result[["table"]], file = "metaplastic_tnbc_nanostringdiff.csv")

###Make a matrix to prepare for normalization 
###that only considers subtype
pheno=pData(NanoStringData)
group=pheno$subtype
design.full=model.matrix(~0+group)
design.full

###perform glm ratio test of mes to average of all other MpBC samples
##note group order, by default, is alphabetical by phenotype label
###so mesenchymal gets +1 designation, while all other 3 control groups add up to -1
###only comparing mesenchymal groups, so TNBC is left out with a '0'
result=glm.LRT(NanoStringData,design.full,contrast=c(1, -1/3, -1/3, -1/3, 0))
head(result$table)
str(result)

###write result out to csv
write.csv(result[["table"]], file = "subtype_mes_vs_metaplastic_nanostringdiff.csv")

###perform glm ratio test of spindle to average of all other MpBC
result=glm.LRT(NanoStringData,design.full,contrast=c(-1/3, -1/3, 1, -1/3, 0))
head(result$table)
str(result)

###write result out to csv
write.csv(result[["table"]], file = "subtype_spindle_vs_metaplastic_nanostringdiff.csv")

###perform glm ratio test of squamous to average of all other MpBC
result=glm.LRT(NanoStringData,design.full,contrast=c(-1/3, -1/3, -1/3, 1, 0))
head(result$table)
str(result)

###write result out to csv
write.csv(result[["table"]], file = "subtype_squamous_vs_metaplastic_nanostringdiff.csv")
