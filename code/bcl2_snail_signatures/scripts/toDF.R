toDF <- function(ngene, wd){
  dirLst <- dir(wd)[grep(paste('_', ngene, '_', sep=""), dir(wd))]
  
  readin_Lst <- list()
  for(i in 1:length(dirLst)){
    readin_Lst[[i]] <- read.csv(paste(wd, dirLst[i], '/pathway_activity_testset.csv', sep=""),
                                row.names="X")
  }
  
  readin_df <- as.data.frame(do.call(cbind, readin_Lst))
  return(readin_df)
}

