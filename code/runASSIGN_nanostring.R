library(ASSIGN)
library(matrixStats)
library(genefilter)
# load input data
indata <- readRDS("data/nanostring_input_data.rds")
geneLists <- readRDS("data/nanostring_panelgeneLists.rds")

runassignGFRN <- function(indata, geneLists, run=c("akt", "bad", "her2", "igf1r",
                                        "krasgv", "raf", "bim", "snail"),
                          optimized_geneList=NULL, use_seed=1234,
                          sigma_sZero=0.05, sigma_sNonZero=0.5,
                          S_zeroPrior=FALSE, iter=100000, burn_in=50000,
                          exclude_common_genes=FALSE, adaptive_S=TRUE, ECM=FALSE) {
  
  #list of anchor genes
  anchorGeneList <- list(akt = "AKT1", bad = "BAD", egfr="EGFR", her2 = "ERBB2",
                         igf1r = "IGF1R", krasgv = "KRAS",
                         raf = "RAF1", bim="BCL2L11", snail=FALSE, becn=FALSE)
  
  #list of corresponding controls for each pathway
  gfpList <- list(akt = "GFP_FOR_AKT", bad = "GFP_18HR", her2 = "GFP_18HR",
                  igf1r = "GFP_18HR", krasgv = "GFP_36HR",
                  raf = "GFP_18HR", bim="GFP_18HR", snail="GFP_18HR")
  
  col_category <- list(akt = "AKT", bad = "BAD", egfr=NULL, her2 = "HER2",
                  igf1r = "IGF1R", krasgv = "KRAS",
                  raf = "RAF1", bim="BIM", snail="SNAIL", becn=NULL)
  

    for (curr_path in run){
      trainingLabel <- list()
      trainingLabel[["control"]] <- list()
      trainingLabel[["control"]][[curr_path]] <- 1:3
      trainingLabel[[curr_path]] <- 4:6
      
      if (!(anchorGeneList[curr_path] %in% rownames(indata))){
        warning(anchorGeneList[curr_path], " not in input data. No anchor gene ",
                "will be used.")
        anchorGeneList[curr_path] <- list(NULL)
      }
      
      excludeGeneList <- NULL
      if (exclude_common_genes){
        excludegenes <- get("excludegenes", envir = environment())
        excludeGeneList <- list()
        excludeGeneList[curr_path] <- list(excludegenes)
      }
      
      if (use_seed){
        set.seed(use_seed)
      }
      
      temp_assay <- indata
      
      if(curr_path %in% c("egfr", "becn")){
        trainingLabel <- NULL
        message("no training data for ", curr_path)
        temp_train <- NULL
        temp_test <- temp_assay[,c(grep(paste0("^","YY"), colnames(indata)), grep("^\\d", colnames(indata)))]
        temp_test <- temp_test[rowSds(temp_test) != 0,]
        message("Test columns: ", paste(colnames(temp_test), collapse = ' '))
        anchorg <- NULL
      } else {
        temp_train <- cbind(temp_assay[,grep(paste0("^",gfpList[[curr_path]]), colnames(indata))],
                            temp_assay[,grep(paste0("^",col_category[[curr_path]]), colnames(indata))])
        if(ncol(temp_train) != 6){
          stop("malformed train ", colnames(temp_train))
        }
        message("Training columns: ", paste(colnames(temp_train), collapse = ' '))
        temp_train <- temp_train[rowSds(temp_train) != 0,]
        temp_test <- temp_assay[,grep(paste0("^","COH"), colnames(indata))]
        temp_test <- temp_test[rowSds(temp_test) != 0,]
        m <- merge_drop(temp_train, temp_test)
        temp_train <- temp_train[rownames(m),]
        temp_test <- temp_test[rownames(m),]
        message("Test columns: ", paste(colnames(temp_test), collapse = ' '))
        anchorg <- anchorGeneList[curr_path]
      }
      tempgene <- list()
      tempgene[curr_path] <- list(geneLists[[curr_path]][geneLists[[curr_path]] %in% rownames(temp_test)])
      message("gene list length: ", length(tempgene[[curr_path]]))
      
      assign.wrapper(trainingData = temp_train,
                     testData = temp_test,
                     anchorGenes = anchorg,
                     excludeGenes = excludeGeneList,
                     trainingLabel = trainingLabel,
                     testLabel = NULL,
                     geneList = tempgene,
                     adaptive_B = TRUE,
                     adaptive_S = adaptive_S,
                     mixture_beta = FALSE,
                     S_zeroPrior = S_zeroPrior,
                     outputDir = paste0(curr_path, "_gene_list_adaptS_",adaptive_S),
                     sigma_sZero = sigma_sZero, sigma_sNonZero = sigma_sNonZero,
                     iter = iter, burn_in = burn_in, ECM = ECM)
    }
}

merge_drop<-function(x,y,by=0,...){
  new_m<-merge(x,y,by=by,...)
  rownames(new_m)<-new_m$Row.names
  return(new_m[,2:length(colnames(new_m))])
}

adapS <- FALSE
setwd("results/")
runassignGFRN(indata = indata, geneLists = geneLists, run="akt", adaptive_S = adapS)
runassignGFRN(indata = indata, geneLists = geneLists, run="bad", adaptive_S = adapS)
runassignGFRN(indata = indata, geneLists = geneLists, run="igf1r", adaptive_S = adapS)
runassignGFRN(indata = indata, geneLists = geneLists, run="her2", adaptive_S = adapS)
runassignGFRN(indata = indata, geneLists = geneLists, run="krasgv", adaptive_S = adapS)
runassignGFRN(indata = indata, geneLists = geneLists, run="raf", adaptive_S = adapS)
runassignGFRN(indata = indata, geneLists = geneLists, run="bim", adaptive_S = adapS)
runassignGFRN(indata = indata, geneLists = geneLists, run="snail", adaptive_S = adapS)
setwd("../")
