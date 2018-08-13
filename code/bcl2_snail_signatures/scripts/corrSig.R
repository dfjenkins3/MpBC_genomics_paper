rm(list=ls())
setwd("/restricted/projectnb/pathsig/yuqingz/emtSig_v2/Incl_adSF/")
source("../toDF.R")
ngene_seq <- c(10, 25, seq(from=0, to=500, by=50)[-1])
geneLst <- c("BCL2L11", "BECN", "SNAI", "Twist")


######## ICBP ########
test_dir <- "/restricted/projectnb/pathsig/20150929_bild_paper_new_ASSIGN/data/"
ICBP_prote <- read.table(paste(test_dir, "correlation_data/proteomics.txt", sep=""))
ICBP_drugres <- read.delim(paste(test_dir, "correlation_data/ICBP_drugs.txt", sep=""), as.is=TRUE)


BCL2L11_prote <- c("cleaved.caspase.7", "Bcl2")
BECN_prote <- c("AMPKp", "mTor") 
SNAI_prote <- c("NOTCH3", "B.Catenin", "stat3", "stat3p705", "E.cadherin")
Twist_prote <- c("FGFR1", "E.cadherin")
ICBP_proteVec <- c(BCL2L11_prote, BECN_prote, SNAI_prote, Twist_prote)
ICBP_proteInd <- c(rep("BCL2L11", length(BCL2L11_prote)),
                   rep("BECN", length(BECN_prote)),
                   rep("SNAI", length(SNAI_prote)),
                   rep("Twist", length(Twist_prote)))

BCL2L11_drugs <- c("Paclitaxel", "Ixabepilone", "Vinorelbine", "Docetaxel")
BECN_drugs <- c("MLN4924", "Rapamycin", "Everolimus", "Temsirolimus")
SNAI_drugs <- c("Disulfiram")
ICBP_drugVec <- c(BCL2L11_drugs, BECN_drugs, SNAI_drugs)
ICBP_drugInd <- c(rep("BCL2L11", length(BCL2L11_drugs)),
                  rep("BECN", length(BECN_drugs)),
                  rep("SNAI", length(SNAI_drugs)))


ICBP_prote_cormat <- matrix(0, nrow=length(ngene_seq), ncol=length(ICBP_proteVec),
                            dimnames=list(ngene_seq, ICBP_proteVec))
ICBP_drug_cormat <- matrix(0, nrow=length(ngene_seq), ncol=length(ICBP_drugVec),
                           dimnames=list(ngene_seq, ICBP_drugVec))

for(ngene in ngene_seq){
  pred_scores <- toDF(ngene, "ICBP/")
  
  ## Compute correlation: Proteins
  # match samples
  samnames_prote <- rownames(ICBP_prote)
  samnames_prote[which(samnames_prote=="SUM1315MO2")] = "SUM1315"
  samnames_prote[which(samnames_prote=="MDAMB157")] = "MB157"
  samnames_pred <- rownames(pred_scores)
  overlap_samnames <- intersect(samnames_prote, samnames_pred)
  ICBP_prote_new <- ICBP_prote[match(overlap_samnames, samnames_prote), ]
  pred_scores_new <- pred_scores[match(overlap_samnames, samnames_pred), ]
  
  # compute correlation
  for(k in 1:length(ICBP_proteVec)){
    ICBP_prote_cormat[as.character(ngene), k] <- cor.test(pred_scores_new[, ICBP_proteInd[k]], 
                                                          ICBP_prote_new[, ICBP_proteVec[k]],
                                                          use="pairwise", method="spearman")$estimate
  }
  
  
  ## Compute correlation: drugs
  # match samples
  samnames_drug <- ICBP_drugres$Cell.line
  
  samnames_pred <- rownames(pred_scores)
  samnames_pred[grep("^X",samnames_pred)] <- gsub("X", "", samnames_pred[grep("^X",samnames_pred)])
  samnames_pred[which(samnames_pred=="MDAMB175")] <- "MDAMB175VII"
  samnames_pred[which(samnames_pred=="SUM1315")] <- "SUM1315MO2"
  samnames_pred[which(samnames_pred=="T47D.Kbluc")] <- "T47D_KBluc"
  
  overlap_samnames <- intersect(samnames_drug, samnames_pred)
  ICBP_drug_new <- ICBP_drugres[match(overlap_samnames, samnames_drug), ]
  pred_scores_new <- pred_scores[match(overlap_samnames, samnames_pred), ]
  
  # compute correlation
  for(k in 1:length(ICBP_drugVec)){
    vec1 <- pred_scores_new[, ICBP_drugInd[k]]
    vec2 <- ICBP_drug_new[, ICBP_drugVec[k]]
    keep <- (!is.na(vec1)) & (!is.na(vec2))
    ICBP_drug_cormat[as.character(ngene), k] <- cor.test(vec1[keep], vec2[keep],
                                                         use="pairwise", method="spearman")$estimate
  }
}

write.csv(ICBP_prote_cormat[, which(ICBP_proteInd=="BCL2L11")], file="../corr_mat_IncladSF/BCL2L11_ICBP_prote.csv")
write.csv(ICBP_prote_cormat[, which(ICBP_proteInd=="BECN")], file="../corr_mat_IncladSF/BECN_ICBP_prote.csv")
write.csv(ICBP_prote_cormat[, which(ICBP_proteInd=="SNAI")], file="../corr_mat_IncladSF/SNAI_ICBP_prote.csv")
write.csv(ICBP_prote_cormat[, which(ICBP_proteInd=="Twist")], file="../corr_mat_IncladSF/Twist_ICBP_prote.csv")

write.csv(ICBP_drug_cormat[, which(ICBP_drugInd=="BCL2L11")], file="../corr_mat_IncladSF/BCL2L11_ICBP_drug.csv")
write.csv(ICBP_drug_cormat[, which(ICBP_drugInd=="BECN")], file="../corr_mat_IncladSF/BECN_ICBP_drug.csv")
write.csv(ICBP_drug_cormat[, which(ICBP_drugInd=="SNAI")], file="../corr_mat_IncladSF/SNAI_ICBP_drug.csv")



######## TCGA ########
rm(list=ls())
source("../toDF.R")
ngene_seq <- c(10, 25, seq(from=0, to=500, by=50)[-1])

test_dir <- "/restricted/projectnb/pathsig/20150929_bild_paper_new_ASSIGN/data/"
TCGA_prote <- read.csv(paste(test_dir, "Datasets/TCGA-BRCA-RBN.csv", sep=""), as.is=TRUE)

BCL2L11_prote <- c("Bim", "Bid", "Bax", "Bak", "Caspase.7_cleavedD198", "Bcl.2", "Bcl.xL")
BECN_prote <- c("Beclin", "mTOR", "AMPK_alpha", "AMPK_pT172") 
SNAI_prote <- c("Smad3", "Smad4", "Notch1", "beta.Catenin", "STAT3_pY705", "Fibronectin", "E.Cadherin", "N.Cadherin")
Twist_prote <- c("TAZ", "Claudin.7", "E.Cadherin", "N.Cadherin")
TCGA_proteVec <- c(BCL2L11_prote, BECN_prote, SNAI_prote, Twist_prote)
TCGA_proteInd <- c(rep("BCL2L11", length(BCL2L11_prote)),
                   rep("BECN", length(BECN_prote)),
                   rep("SNAI", length(SNAI_prote)),
                   rep("Twist", length(Twist_prote)))

TCGA_prote_cormat <- matrix(0, nrow=length(ngene_seq), ncol=length(TCGA_proteVec),
                            dimnames=list(ngene_seq, TCGA_proteVec))

for(ngene in ngene_seq){
  pred_scores <- toDF(ngene, "TCGA/")
  
  ## Match samples
  samnames_prote <- TCGA_prote[, "TCGA_patient_barcode"]
  samnames_prote <- gsub('-','.',samnames_prote)
  
  samnames_pred <- rownames(pred_scores)
  for(i in 1:length(samnames_pred)){
    samnames_pred[i] <- substr(samnames_pred[i], 1, 12)
  }
  
  overlap_samnames <- intersect(samnames_prote, samnames_pred)
  TCGA_prote_new <- TCGA_prote[match(overlap_samnames, samnames_prote), ]
  pred_scores_new <- pred_scores[match(overlap_samnames, samnames_pred), ]
  
  
  ## Compute correlation
  for(k in 1:length(TCGA_proteVec)){
    vec1 <- pred_scores_new[, TCGA_proteInd[k]]
    vec2 <- TCGA_prote_new[, TCGA_proteVec[k]]
    keep <- (!is.na(vec1)) & (!is.na(vec2))
    #print(sum(!keep))
    TCGA_prote_cormat[as.character(ngene), k] <- cor.test(vec1[keep], vec2[keep],
                                                          use="pairwise", method="spearman")$estimate
  }
  
}    

write.csv(TCGA_prote_cormat[, which(TCGA_proteInd=="BCL2L11")], file="../corr_mat_IncladSF/BCL2L11_TCGA_prote.csv")
write.csv(TCGA_prote_cormat[, which(TCGA_proteInd=="BECN")], file="../corr_mat_IncladSF/BECN_TCGA_prote.csv")
write.csv(TCGA_prote_cormat[, which(TCGA_proteInd=="SNAI")], file="../corr_mat_IncladSF/SNAI_TCGA_prote.csv")
write.csv(TCGA_prote_cormat[, which(TCGA_proteInd=="Twist")], file="../corr_mat_IncladSF/Twist_TCGA_prote.csv")
