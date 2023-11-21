pgg <- function(
    phenoDTfile= NULL,
    analysisId=NULL,
    trait=NULL, # per trait
    environment=NULL,
    proportion=10,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATES THE PREDICTED GENETIC GAIN FOR TRAITS
  ## IS USED IN THE BANAL APP UNDER THE METRICS MODULES
  pggAnalysisId <- as.numeric(Sys.time())
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the ID of the analysis to use as input", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  
  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions 
  mydata <- mydata[which(mydata$analysisId %in% analysisId),]
  if(nrow(mydata)==0){stop("No match for this analysisId. Please correct.", call. = FALSE)}
  if(is.null(environment)){environment <- na.omit(unique(mydata$environment))}
  if(is.null(phenoDTfile$data$pedigree) || (nrow(phenoDTfile$data$pedigree) == 0 ) ){stop("yearOfOrigin column was not matched in your original file. Please correct.", call. = FALSE)}
  yearsToUse <- as.character(unique(phenoDTfile$data$pedigree$yearOfOrigin))
  mydata <- merge(mydata, phenoDTfile$data$pedigree[,c("designation","yearOfOrigin")], by="designation", all.x=TRUE )
  mydata <- mydata[which(!is.na(mydata$yearOfOrigin)),]
  yearsTesting <- unique(phenoDTfile$data$pheno[,c("designation","year")])
  yearsTesting <- yearsTesting[which(!duplicated(yearsTesting$designation)),]
  mydata <- merge(mydata, yearsTesting, by="designation", all.x=TRUE )
  mydata <- mydata[which(!is.na(mydata$year)),]
  ############################
  ## gg analysis
  p <- proportion/100
  i <- dnorm(qnorm(1 - p))/p
  # counter=1
  for(iTrait in trait){ # iTrait = trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    uEnvironments <- unique(mydata$environment)
    for(uE in uEnvironments){ # uE <- uEnvironments[1]
      # subset data
      mydataSub <- droplevels(mydata[which((mydata$trait == iTrait) & (mydata$environment %in% uE)  ),])
      # calculate parameters
      rels <- mydataSub$rel;
      badrels <- which(rels < 0); if(length(badrels) > 0){ rels[badrels] <- 1e-6}
      r <- ifelse(length(na.omit(rels)) > 0, mean(sqrt(na.omit(rels)), na.rm=TRUE), 1e-6)
      sigma<- sd(mydataSub$predictedValue, na.rm = TRUE)
      mu<- mean(mydataSub$predictedValue, na.rm = TRUE)
      mydataSubSorted <- mydataSub[with(mydataSub, order(-predictedValue)), ]
      mydataSubSortedSel <- mydataSubSorted[1:round(nrow(mydataSubSorted) * p),]
      age <- mean(mydataSubSortedSel$year, na.rm=TRUE) - mean(mydataSubSortedSel$yearOfOrigin, na.rm=TRUE)
      R <- r * sigma * i
      ggAge =  R/age
      ##
      phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                   data.frame(module="pgg",analysisId=pggAnalysisId, trait= iTrait, environment=uE, 
                                              parameter=c("r","sigmaG","meanG", "cycleLength","i","R","PGG"), method=c("sqrt(r2)","sd(BLUP)","sum(x)/n","yearTest-yearOrigin","dnorm(qnorm(1 - p))/p","r*sigma*i","R/cycleLength"), 
                                              value=c(r,sigma, mu, age, i, R, ggAge), 
                                              stdError=0
                                   )
      )
      currentModeling <- data.frame(module="pgg", analysisId=pggAnalysisId,trait=iTrait, environment=uE, 
                                    parameter=c("proportion","verbose"), value=c(proportion, verbose))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
    }
    
  }
  #########################################
  # update databases
  phenoDTfile$status <- rbind( phenoDTfile$status, data.frame(module="pgg", analysisId=pggAnalysisId))
  ## add which data was used as input
  modeling <- data.frame(module="pgg",  analysisId=pggAnalysisId, trait=c("inputObject"), environment="general", 
                         parameter= c("analysisId"), value= c(analysisId ))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)#paste("pgg done:",id))
}
