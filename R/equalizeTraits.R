equalizeTraits <- function(object, traits, newName=NULL){
  mydata <- object$data$pheno
  newValue <- apply(mydata[,traits, drop=FALSE],1, function(x){mean(x, na.rm=TRUE)})
  if(is.null(newName) ){
    newName <- paste(traits,collapse = "_")
  }
  if( newName == "" ){
    newName <- paste(traits,collapse = "_")
  }
  mydata[,newName] <- newValue
  object$data$pheno <- mydata
  object$metadata$pheno <- rbind( object$metadata$pheno, data.frame(parameter="trait", value=newName))
  return(object)
}
