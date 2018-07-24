library(ASSIGN)
library(matrixStats)
library(genefilter)

# load input data
geneLists <- readRDS("data/nanostring_geneLists.rds")
anastassiou <- read.table("data/Anastassiou.txt", header=T, stringsAsFactors = F)
boquest_stemcell_up <- read.table("data/Boquest_stemcell_up.txt", header=T, stringsAsFactors = F)
boquest_stemcell_dn <- read.table("data/Boquest_stemcell_dn.txt", header=T, stringsAsFactors = F)
phong_tnf_up <- read.table("data/Phong_TNF_up.txt", header=T, stringsAsFactors = F)
indataall <- readRDS("data/nanostring_input_data.rds")
temp_test_all <- indataall[,grep(paste0("^","COH"), colnames(indataall))]
temp_test_all <- temp_test_all[rowSds(temp_test_all) != 0,]

#anastassiou signature, no training data, all up
anastassiou_direction <- rep(1,25)
assign.wrapper(trainingData = NULL,
               testData = temp_test_all,
               anchorGenes = NULL,
               excludeGenes = NULL,
               trainingLabel = NULL,
               testLabel = NULL,
               geneList = geneLists[1],
               adaptive_B = TRUE,
               adaptive_S = TRUE,
               mixture_beta = FALSE,
               S_zeroPrior = FALSE,
               outputDir = paste0("results/anastassiou_gene_list_adaptS_TRUE_ALL"),
               sigma_sZero = 0.05, sigma_sNonZero = 0.5,
               iter = 100000, burn_in = 50000, ECM = FALSE,
               override_S_matrix = anastassiou_direction)

#Boquest Stem Cell signature, no training data, 42 up, 12 down
tester <- assign.preprocess(trainingData = NULL,
                            trainingLabel = NULL,
                            testData = temp_test_all,
                            geneList = geneLists[2])
boquest_stemcell_direction <- ifelse(rownames(tester$testData_sub) %in% boquest_stemcell_up$Boquest_stemcell_up, 1, -1)
assign.wrapper(trainingData = NULL,
               testData = temp_test_all,
               anchorGenes = NULL,
               excludeGenes = NULL,
               trainingLabel = NULL,
               testLabel = NULL,
               geneList = geneLists[2],
               adaptive_B = TRUE,
               adaptive_S = TRUE,
               mixture_beta = FALSE,
               S_zeroPrior = FALSE,
               outputDir = paste0("results/boquest_stemcell_gene_list_adaptS_TRUE_ALL"),
               sigma_sZero = 0.05, sigma_sNonZero = 0.5,
               iter = 100000, burn_in = 50000, ECM = FALSE,
               override_S_matrix = boquest_stemcell_direction)

#Phong TNF signature, no training data, all genes expected to increase
phong_tnf_direction <- rep(1,27)
assign.wrapper(trainingData = NULL,
               testData = temp_test_all,
               anchorGenes = NULL,
               excludeGenes = NULL,
               trainingLabel = NULL,
               testLabel = NULL,
               geneList = geneLists[3],
               adaptive_B = TRUE,
               adaptive_S = TRUE,
               mixture_beta = FALSE,
               S_zeroPrior = FALSE,
               outputDir = paste0("results/phong_tnf_gene_list_adaptS_TRUE_ALL"),
               sigma_sZero = 0.05, sigma_sNonZero = 0.5,
               iter = 100000, burn_in = 50000, ECM = FALSE,
               override_S_matrix = phong_tnf_direction)
