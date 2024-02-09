
singleCrossMat <- function( # single cross matrix function
  object= NULL,
  hybridBatch=1000,
  allHybrids=FALSE,
  verbose=TRUE
){
  analysisId <- as.numeric(Sys.time())
  ############################
  # loading the dataset
  if (is.null(object)) stop("No input file specified.")
  
  ## extract marker matrices and reference alleles
  M1 <- object$data$geno
  M2 <- object$data$geno
  
  ############################
  ## first get all possible hybrids
  if(allHybrids){ # reduce to only hybrids present in the designation column in the pedigree table
    res1 <- sommer::build.HMM(M1, M2, 
                              return.combos.only = TRUE)
  }else{
    ped <- object$data$pedigree
    metaPed <- object$metadata$pedigree
    colnames(ped) <- cgiarBase::replaceValues(colnames(ped), Search = metaPed$value, Replace = metaPed$parameter )
    cross <- unique(ped[,c("mother","father","designation")]); colnames(cross) <- c("Var1","Var2","hybrid")
    possible <- apply(cross,1, function(x){
      ifelse(x[1]%in%rownames(M1),1,0) + ifelse(x[2]%in%rownames(M1),1,0)
    })
    possible <- which(possible == 2)
    res1 <- list(data.used=cross[possible,]) 
  }
  
  if(nrow(res1$data.used)>0){ # if there is hybrids to build
    ## build the marker matrix for batches of 1000 hybrids
    batches <- sort(rep(1:1000,min(c(nrow(res1$data.used),hybridBatch))))
    res1$data.used$batch <- batches[1:nrow(res1$data.used)]
    data.usedBatches <- split(res1$data.used, res1$data.used$batch)
    # start the loop
    for(i in 1:length(data.usedBatches)){
      prov <- sommer::build.HMM(M1, M2,
                                custom.hyb = data.usedBatches[[i]], 
                                separator=separator
      )
      if(i == 1){
        M <- prov$HMM.add
      }else{
        M <- rbind(M,prov$HMM.add)
      }
    }
  }else{
    M <- matrix(NA, nrow=0, ncol=ncol(M1)); colnames(M) <- colnames(M1)
  }
  object$data$geno <- rbind(object$data$geno,M)
  # other tables
  object$status <- rbind( object$status, data.frame(module="scm", analysisId=analysisId))
  ## add which data was used as input
  modeling <- data.frame(module="scm",  analysisId=analysisId, trait=c("none"), environment="general",
                         parameter= c("allHybrids"), value= c(allHybrids ))
  if(is.null(object$modeling)){
    object$modeling <-  modeling
  }else{
    object$modeling <- rbind(object$modeling, modeling[, colnames(object$modeling)])
  }
  
  ##
  object$metrics <- rbind(object$metrics,
                               data.frame(module="scm",analysisId=analysisId, trait= "none", environment="across",
                                          parameter=c("nHybrids"), method=c("(Md + Mf) * (1/2)"),
                                          value=c( nrow(M)  ),
                                          stdError=c(NA)
                               )
  )
  return(object)
}
