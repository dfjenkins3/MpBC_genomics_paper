# grid search through correlation results
# run ASSIGN and protein correlation as described in
# https://doi.org/10.1186/s13073-017-0429-x to
# generate results for all the pathway lengths listed below
library(ASSIGN)

results <- expand.grid(akt=c(5,10,15,20),
                       bad=c(5,10,15,20,25,50,75,100,125,150,175,200,225,250),
                       egfr=c(5,10,15,20,25,50),
                       her2=c(5,10),
                       igf1r=c(5,10,15,20,25,50,75,100),
                       krasgv=c(5,10,15,20,25,50,75,100,125,150,175,200),
                       raf=c(5,10,15,20,25,50,75,100,125,150,175,200,225,250,275,300,350))
results$sum_withoverlap <- rowSums(results)
data("gfrn_geneList")

num_uniquegenes <- function(akt, bad, egfr, her2, igf1r, krasgv, raf, genelist){
  length(unique(c(genelist[["akt_up"]][1:floor(akt/2)],
                  genelist[["akt_down"]][1:ceiling(akt/2)],
                  genelist[["bad_up"]][1:floor(bad/2)],
                  genelist[["bad_down"]][1:ceiling(bad/2)],
                  genelist[["egfr_up"]][1:floor(egfr/2)],
                  genelist[["egfr_down"]][1:ceiling(egfr/2)],
                  genelist[["her2_up"]][1:floor(her2/2)],
                  genelist[["her2_down"]][1:ceiling(her2/2)],
                  genelist[["igf1r_up"]][1:floor(igf1r/2)],
                  genelist[["igf1r_down"]][1:ceiling(igf1r/2)],
                  genelist[["krasgv_up"]][1:floor(krasgv/2)],
                  genelist[["krasgv_down"]][1:ceiling(krasgv/2)],
                  genelist[["raf_up"]][1:floor(raf/2)],
                  genelist[["raf_down"]][1:ceiling(raf/2)])))
}

results$sum_nooverlap <- apply(results, 1,
                               function(x) num_uniquegenes(x[1], x[2], x[3], x[4], x[5], x[6], x[7], gfrn_geneList))
saveRDS(results, file="results/parametermat.rds")

icbp_protein <- read.table("data/icbp/cor_prot.txt", sep="\t", header=T, row.names=1)
icbp_drug <- read.table("data/icbp/cor_drug_mat.txt", sep="\t", header=T, row.names=1)
tcga_protein <- read.table("data/tcga/cor.txt", sep="\t", row.names=1, header=T)

get_cor_value <- function(corfile, pathway, genelist, protein){
  corfile[paste(paste(pathway,as.character(genelist),"gene","list",sep="_"), "adap_adap_single", pathway, sep="/"),protein]
}

#akt
message("akt")
##icbp Akt
results$icbp_akt_Akt <- apply(results, 1, function(x) get_cor_value(icbp_protein,"akt",x[1], "Akt"))
##icbp PDK1
system.time(results$icbp_akt_PDK1 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"akt",x[1], "PDK1")))
##icbp PDK1p241
results$icbp_akt_PDK1p241 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"akt",x[1], "PDK1p241"))
#bad
message("bad")
##icbp negAkt
results$icbp_bad_negAkt <- -1 * apply(results, 1, function(x) get_cor_value(icbp_protein,"bad",x[2], "Akt"))
##icbp negPDK1
results$icbp_bad_negPDK1 <- -1 * apply(results, 1, function(x) get_cor_value(icbp_protein,"bad",x[2], "PDK1"))
##icbp negPDK1p241
results$icbp_bad_negPDK1p241 <- -1 * apply(results, 1, function(x) get_cor_value(icbp_protein,"bad",x[2], "PDK1p241"))
#egfr
message("egfr")
##icbp EGFR
results$icbp_egfr_EGFR <- apply(results, 1, function(x) get_cor_value(icbp_protein,"egfr",x[3], "EGFR"))
##icbp EGFRp1068
results$icbp_egfr_EGFRp1068 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"egfr",x[3], "EGFRp1068"))
#her2
message("her2")
##icbp HER2
results$icbp_her2_HER2 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"her2",x[4], "HER2"))
##icbp HER2p1248
results$icbp_her2_HER2p1248 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"her2",x[4], "HER2p1248"))
#igf1r
message("igf1r")
##icbp IGFR1
results$icbp_igf1r_IGFR1 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"igf1r",x[5], "IGFR1"))
##icbp PDK1
results$icbp_igf1r_PDK1 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"igf1r",x[5], "PDK1"))
##icbp PDK1p241
results$icbp_igf1r_PDK1p241 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"igf1r",x[5], "PDK1p241"))
#krasgv
message("krasgv")
##icbp EGFR
results$icbp_krasgv_EGFR <- apply(results, 1, function(x) get_cor_value(icbp_protein,"krasgv",x[6], "EGFR"))
##icbp EGFRp1068
results$icbp_krasgv_EGFRp1068 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"krasgv",x[6], "EGFRp1068"))
#raf
message("raf")
##icbp MEK1
results$icbp_raf_MEK1 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"raf",x[7], "MEK1"))
results[,25] <- apply(results, 1, function(x) get_cor_value(icbp_protein,"raf",x[7], "MEK1"))
colnames(results)[25] <- "icbp_raf_MEK1"
##icbp PKCalphap657
results$icbp_raf_PKCalphap657 <- apply(results, 1, function(x) get_cor_value(icbp_protein,"raf",x[7], "PKCalphap657"))
##icbp PKCalpha
results$icbp_raf_PKCalpha <- apply(results, 1, function(x) get_cor_value(icbp_protein,"raf",x[7], "PKCalpha"))

save.image(file="results/session_done.RData")
