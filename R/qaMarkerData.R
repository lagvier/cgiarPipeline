qaMarker <- function(
    markerDTfile = NULL, # markers are assumed to come centered
    propNaUpperThreshForMarker=.3,
    propNaUpperThreshForInds=.3,
    maf=0.05, ploidy=2,
    imputationMethod="median"
){
  ## THIS FUNCTION IMPUTES A MARKER MATRIX
  ## IS USED IN THE BCLEAN APP UNDER THE QA/QC MODULES
  qamAnalysisId <- as.numeric(Sys.time())
  modifications <- as.data.frame(matrix(nrow=0, ncol=6))
  colnames(modifications) <- c("module" ,     "analysisId"  ,"reason"  ,     "row" ,"col"  , "value" )
  ############################
  # loading the dataset
  M <- markerDTfile$data$geno
  ## identify which individuals should be removed according to the threshold
  propNaIndividual <- apply(M,1,function(x){(length(which(is.na(x)))/length(x))})
  badIndividual <- which(propNaIndividual > propNaUpperThreshForInds)
  if(length(badIndividual) > 0){
    modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="%missing", row=badIndividual, col=NA, value=NA) )
  }
  ## identify which markers should be removed according to the threshold
  propNaMarker <- apply(M,2,function(x){(length(which(is.na(x)))/length(x))})
  badMarker <- which(propNaMarker > propNaUpperThreshForMarker)
  if(length(badMarker) > 0){
    modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="%missing", row=NA, col=badMarker, value=NA) )
  }
  ## identify bad MAF
  MAF <- apply(M+1, 2, function(x) {
    AF <- mean(x, na.rm = T)/ploidy
    MAF <- ifelse(AF > 0.5, 1 - AF, AF)
  })
  badMarker2 <- which(MAF < maf)
  if(length(badMarker2) > 0){
    modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="MAF", row=NA, col=badMarker2, value=NA) )
  }
  # imputation track
  toImpute <- which(is.na(M), arr.ind = TRUE)
  remove1 <- which(toImpute[,1] %in% badIndividual)
  if(length(remove1) > 0){toImpute <- toImpute[-remove1,]}
  remove2 <- which(toImpute[,2] %in% badMarker)
  if(length(remove2) > 0){toImpute <- toImpute[-remove2,]}
  remove3 <- which(toImpute[,2] %in% badMarker2)
  if(length(remove3) > 0){toImpute <- toImpute[-remove3,]}
  # impute
  if(imputationMethod == "median"){
    M <- apply(M,2,sommer::imputev)
    if(nrow(toImpute) > 0){
      modifications <- rbind(modifications, data.frame(module="qaGeno", analysisId=qamAnalysisId, reason="impute", row=toImpute[,1], col=toImpute[,2], value=M[toImpute] ) )
    }
  }else{
    stop("Method not implemented yet.",call. = FALSE)
  }
  markerDTfile$modifications$geno <- rbind(markerDTfile$modifications$geno, modifications)
  #########################################
  ## update the rest of the data tables
  ## write the parameters to the parameter database
  markerDTfile$status <- rbind( markerDTfile$status, data.frame(module="qaGeno", analysisId=qamAnalysisId))

  return(markerDTfile)
}
