# grid search through correlation results TCGA

load("results/session_done.RData")

#akt
message("akt")
##tcga Akt
results$tcga_akt_Akt <- apply(results, 1, function(x) get_cor_value(tcga_protein,"akt",x[1], "Akt"))
##tcga PDK1
results$tcga_akt_PDK1 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"akt",x[1], "PDK1"))
##tcga PDK1p241
results$tcga_akt_PDK1_pS241 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"akt",x[1], "PDK1_pS241"))
#bad
message("bad")
##tcga negAkt
results$tcga_bad_negAkt <- -1 * apply(results, 1, function(x) get_cor_value(tcga_protein,"bad",x[2], "Akt"))
##tcga Bad_pS112
results$tcga_bad_Bad_pS112 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"bad",x[2], "Bad_pS112"))
##tcga negPDK1
results$tcga_bad_negPDK1 <- -1 * apply(results, 1, function(x) get_cor_value(tcga_protein,"bad",x[2], "PDK1"))
##tcga negPDK1p241
results$tcga_bad_negPDK1_pS241 <- -1 * apply(results, 1, function(x) get_cor_value(tcga_protein,"bad",x[2], "PDK1_pS241"))
#egfr
message("egfr")
##tcga EGFR
results$tcga_egfr_EGFR <- apply(results, 1, function(x) get_cor_value(tcga_protein,"egfr",x[3], "EGFR"))
##tcga EGFR_pY1068
results$tcga_egfr_EGFR_pY1068 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"egfr",x[3], "EGFR_pY1068"))
##tcga EGFR_pY1173
results$tcga_egfr_EGFR_pY1173 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"egfr",x[3], "EGFR_pY1173"))
#her2
message("her2")
##tcga HER2
results$tcga_her2_HER2 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"her2",x[4], "HER2"))
##tcga HER2_pY1248
results$tcga_her2_HER2_pY1248 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"her2",x[4], "HER2_pY1248"))
#igf1r
message("igf1r")
##tcga IRS1
results$tcga_igf1r_IRS1 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"igf1r",x[5], "IRS1"))
##tcga PDK1
results$tcga_igf1r_PDK1 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"igf1r",x[5], "PDK1"))
##tcga PDK1p241
results$tcga_igf1r_PDK1_pS241 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"igf1r",x[5], "PDK1_pS241"))
#krasgv
message("krasgv")
##tcga EGFR
results$tcga_krasgv_EGFR <- apply(results, 1, function(x) get_cor_value(tcga_protein,"krasgv",x[6], "EGFR"))
##tcga EGFR_pY1068
results$tcga_krasgv_EGFR_pY1068 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"krasgv",x[6], "EGFR_pY1068"))
##tcga EGFR_pY1173
results$tcga_krasgv_EGFR_pY1173 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"krasgv",x[6], "EGFR_pY1173"))
##tcga MEK1
results$tcga_krasgv_MEK1 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"krasgv",x[6], "MEK1"))
#raf
message("raf")
##tcga MEK1
results$tcga_raf_MEK1 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"raf",x[7], "MEK1"))
##tcga PKC.alpha
results$tcga_raf_PKCalpha <- apply(results, 1, function(x) get_cor_value(tcga_protein,"raf",x[7], "PKC.alpha"))
##tcga PKC.alpha_pS657
results$tcga_raf_PKCalpha_pS657 <- apply(results, 1, function(x) get_cor_value(tcga_protein,"raf",x[7], "PKC.alpha_pS657"))

save.image(file="results/session_done_tcga.RData")