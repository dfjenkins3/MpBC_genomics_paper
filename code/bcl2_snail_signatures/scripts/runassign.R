runassign <- function(indata, run=c("BCL2L11","BECN","SNAI","Twist"), ngene=ngene, 
                      testname=testname, outputpathname="./",
                      optimized_geneList=NULL, exclude_geneList=NULL,
                      use_seed=1234, sigma_sZero=0.05, sigma_sNonZero=0.5,
                      S_zeroPrior=FALSE, iter=100000, burn_in=50000) {
  
  #list of anchor genes
  anchorGeneList <- list(BCL2L11="BCL2L11", BECN="BECN1", SNAI="SNAI1", Twist="TWIST1")
  
  #list of corresponding controls for each pathway
  gfpList <- list(BCL2L11="GFP1", BECN="GFP2", SNAI="GFP1", Twist="GFP1")
  
  for (curr_path in run){
    trainingLabel <- list()
    trainingLabel[['control']] <- list()
    trainingLabel[['control']][[curr_path]] <- 1:
      ncol(indata[[gfpList[[curr_path]]]])
    trainingLabel[[curr_path]] <- (ncol(indata[[gfpList[[curr_path]]]])+1):
      (ncol(indata[[gfpList[[curr_path]]]])+ncol(indata[[curr_path]]))
    
    if(!(anchorGeneList[curr_path] %in% rownames(indata[['test']]))){
      warning(anchorGeneList[curr_path], " not in input data. No anchor gene ",
              "will be used.")
      anchorGeneList[curr_path] <- list(NULL)
    }
    
    if(use_seed){
      set.seed(use_seed)
    }
    
    assign.wrapper(trainingData=cbind(indata[[gfpList[[curr_path]]]],
                                      indata[[curr_path]]),
                   testData=indata[['test']],
                   anchorGenes=anchorGeneList[curr_path],
                   excludeGenes=exclude_geneList[curr_path],
                   trainingLabel=trainingLabel,
                   geneList=NULL,
                   n_sigGene=rep(ngene, 1),
                   adaptive_B=TRUE,
                   adaptive_S=TRUE,
                   mixture_beta=FALSE,
                   S_zeroPrior=S_zeroPrior,
                   outputDir=paste(outputpathname, testname, "/", curr_path,"_", ngene, "_res", sep=""),
                   sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                   iter=iter, burn_in=burn_in)
  }
}