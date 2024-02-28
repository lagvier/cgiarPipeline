equalizeTraits <- function(object, traits, newName=NULL){
  mydata <- object$data$pheno
  newValue <- apply(mydata[,traits, drop=FALSE],1, function(x){mean(x, na.rm=TRUE)})
  if(is.null(newName)){
    mydata[,paste(traits,collapse = "_")] <- newValue
  }else{
    mydata[,newName] <- newValue
  }
  object$data$pheno <- mydata
  return(object)
}
