metLMM <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    analysisId=NULL,
    fixedTerm= c("environment"),  randomTerm=c("designation"),  residualBy=NULL,
    interactionsWithGeno=NULL, envsToInclude=NULL,
    trait= NULL, traitFamily=NULL, useWeights=TRUE,
    heritLB= 0.15,  heritUB= 0.95,
    modelType="blup", # either "grm", "nrm", or both
    deregress=FALSE,  nPC=0,
    maxIters=50, batchSizeToPredict=500, tolParInv=1e-4,
    verbose=TRUE
){
  ## THIS FUNCTION PERFORMS A MULT TRIAL ANALYSIS USING LMM SOLVER
  ## IS USED IN THE BANAL APP UNDER THE GENETIC EVALUATION MODULES
  mtaAnalysisId <- as.numeric(Sys.time())
  if(is.null(phenoDTfile)){stop("Please provide the phenotype file", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the analysisId to be analyzed", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(traitFamily)){traitFamily <- rep("quasi(link = 'identity', variance = 'constant')", length(trait))}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  if(length(fixedTerm) == 0){stop("Please provide your fixed effect terms.", call. = FALSE)}
  if(length(randomTerm) == 0){stop("Please provide your randomTerm effect terms.", call. = FALSE)}
  if(modelType %in% c("gblup","ssgblup","rrblup") & is.null(phenoDTfile$data$geno)){stop("Please include marker information in your data structure to fit this model type.", call. = FALSE)}
  if(is.null(phenoDTfile$metadata$weather)){
    provMet <- as.data.frame(matrix(nrow=0, ncol=3))
    colnames(provMet) <- c("environment", "parameter" ,  "value")
    phenoDTfile$metadata$weather <- provMet
  }
  names(traitFamily) <- trait
  heritLB <- rep(heritLB,length(trait))
  heritUB <- rep(heritUB,length(trait))
  traitOrig <- trait
  common <- intersect(fixedTerm,randomTerm)
  fixedTerm <- setdiff(fixedTerm,common)
  if("designation" %in% randomTerm){}else{interactionsWithGeno <- NULL} # don't allow interactions if genotype is not random from start
  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions #
  if (nrow(mydata) < 2) stop("Not enough data is available to perform a multi trial analysis. Please perform an STA before trying to do an MET.", call. = FALSE)
  mydata <- mydata[which(mydata$analysisId %in% analysisId),]
  if(nrow(mydata)==0){stop("No match for this analysisId. Please correct.", call. = FALSE)}
  if( length(setdiff(setdiff(fixedTerm,"1"),colnames(mydata))) > 0 ){stop(paste("column(s):", paste(setdiff(setdiff(fixedTerm,"1"),colnames(mydata)), collapse = ","),"couldn't be found."), call. = FALSE)}
  if( length(setdiff(setdiff(randomTerm,"1"),colnames(mydata))) > 0 ){stop(paste("column(s):", paste(setdiff(setdiff(randomTerm,"1"),colnames(mydata)), collapse = ","),"couldn't be found."), call. = FALSE)}
  metrics <- phenoDTfile$metrics 
  metrics <- metrics[which(metrics$analysisId %in% analysisId),]
  #############################
  # loading the additional matrix needed depending on the model
  if(modelType == "blup"){
    surrogate <- "TGV"
  }else if(modelType == "pblup"){
    surrogate <- "EBV"
  }else if(modelType %in% c("gblup","rrblup")){
    surrogate <- "GEBV"
  }else if(modelType == "ssgblup"){
    surrogate <- "ssGEBV"
  }
  # if user didn't provide a table for which environments should be included make it
  if(is.null(envsToInclude)){
    envsToInclude=  as.data.frame( do.call( rbind, list (with(mydata, table(environment,trait)) ) ) )
    bad <- which(envsToInclude <= 1, arr.ind = TRUE)
    if(nrow(bad) > 0){envsToInclude[bad] = 0}
    envsToInclude[which(envsToInclude > 1, arr.ind = TRUE)] = 1
  }
  # check if traits specified are actualy available and remove the ones that are not present
  utraits <- unique(mydata$trait)
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% utraits){
      if(verbose){ cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))}
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  traitToRemove <- character()
  if(!is.null(envsToInclude)){
    for(k in 1:length(trait)){
      if( length(which(envsToInclude[,trait[k]] > 0)) == 0 ){
        if(verbose){ cat(paste0("'", trait[k], "' has not envs selected. It will be removed from trait list \n"))}
        traitToRemove <- c(traitToRemove,trait[k])
      }
    }
    trait <- setdiff(trait,traitToRemove)
  }
  if(length(trait)==0){stop("None of the traits specified are available. Please double check", call. = FALSE)}
  heritLB <- heritLB[which(traitOrig %in% trait)]
  heritUB <- heritUB[which(traitOrig %in% trait)]
  ############################
  # remove outliers from each environment
  envs <- unique(mydata$environment)
  mydata$rowindex <- 1:nrow(mydata)
  ##############################
  ## met analysis
  predictionsList <- list(); counter=counter2=1
  traitToRemove <- character()
  for(iTrait in trait){ # # iTrait = trait[1]  iTrait="value"
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    # subset data
    mydataSub <- droplevels(mydata[which(mydata$trait == iTrait),])
    # if(!is.null(envsToInclude)){
    #   envsToIncludeTrait = rownames(envsToInclude)[which(envsToInclude[,iTrait] > 0)]
    #   mydataSub <- mydataSub[which(mydataSub$environment %in% envsToIncludeTrait),]
    # }
    # remove bad environment
    pipeline_metricsSub <- metrics[which(metrics$trait == iTrait & metrics$parameter %in% c("plotH2","H2","meanR2")),]
    goodFields <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value > heritLB[counter2]) & (pipeline_metricsSub$value < heritUB[counter2])),"environment"])
    mydataSub <- mydataSub[which(mydataSub$environment %in% goodFields),]
    if(verbose){print(paste("Fields included:",paste(goodFields,collapse = ",")))}
    ## remove records without marker data if marker effects
    if(modelType == "rrblup"){ # if rrBLUP model we remove records without markers
      mydataSub <- mydataSub[which(mydataSub$designation %in% rownames(phenoDTfile$data$geno)),]
      Mtrait <- phenoDTfile$data$geno[which(rownames(phenoDTfile$data$geno) %in% mydataSub$designation),]
    }
    LGrp <- list();   groupTrait <- NULL
    ## next step
    if(nrow(mydataSub) == 0){
      print(paste("There was no predictions to work with in trait",iTrait,". Please look at your H2 boundaries. You may be discarding all envs."))
      traitToRemove <- c(traitToRemove,iTrait)
    }else{
      mydataSub$designation <- as.factor(mydataSub$designation) # move to factor
      mydataSub$environment <- as.factor(mydataSub$environment) # move to factor
      mydataSub$pipeline <- as.factor(mydataSub$pipeline) # move to factor
      mydataSub <- mydataSub[which(mydataSub$designation != ""),] # remove blank designations
      ## build the environmental index
      ei <- aggregate(predictedValue~environment, data=mydataSub,FUN=mean); colnames(ei)[2] <- "envIndex0"
      ei <- ei[with(ei, order(envIndex0)), ]
      ei$envIndex <- ei$envIndex0 - mean(ei$envIndex0)
      colnames(ei) <- cgiarBase::replaceValues(colnames(ei), Search = "envIndex0", Replace = "value")
      ei$parameter <- paste0("envIndex_",iTrait)
      # update the weather metadata
      phenoDTfile$metadata$weather <- rbind(phenoDTfile$metadata$weather,ei[,colnames(phenoDTfile$metadata$weather)])
      toKeep <- rownames(unique(phenoDTfile$metadata$weather[,c("environment","parameter")]))
      phenoDTfile$metadata$weather <- phenoDTfile$metadata$weather[toKeep,]
      ## add metadata from environment(weather)
      if(!is.null(phenoDTfile$metadata$weather)){
        metas <- phenoDTfile$metadata$weather; 
        metas <- reshape(metas, direction = "wide", idvar = "environment",
                         timevar = "parameter", v.names = "value", sep= "_")
        colnames(metas) <- gsub("value_","", colnames(metas))
        metasClass <- unlist(lapply(metas,class))
        numericMetas <- names(metasClass)[which(metasClass %in% c("integer","numeric"))]
        # center variables
        for(iMeta in numericMetas){
          metas[,iMeta] <- metas[,iMeta] - mean(metas[,iMeta], na.rm=TRUE)
        }
        metas <- metas[which(metas$environment %in% goodFields),]
        colnames(metas) <- cgiarBase::replaceValues(colnames(metas), Search = paste0("envIndex_", iTrait), Replace = "envIndex")
        mydataSub <- merge(mydataSub,metas, by="environment", all.x = TRUE)
      }
      ## define the interactions to be used
      if(!is.null(interactionsWithGeno)){
        interactionsWithGenoTrait <- interactionsWithGeno
        interactionsWithGenoToRemove <- character()
        for(iInter in 1:length(interactionsWithGenoTrait)){
          if(interactionsWithGenoTrait[iInter] %in% colnames(mydataSub)){ # if trait is even there in dataset
            checkInters <- length(unique(mydataSub[,interactionsWithGenoTrait[iInter]]))
            if (checkInters < 2){ # there needs to be at least more than one level
              interactionsWithGenoToRemove <- c(interactionsWithGenoToRemove,interactionsWithGenoTrait[iInter])
            }
          }else{ # remove straight away
            interactionsWithGenoToRemove <- c(interactionsWithGenoToRemove,interactionsWithGenoTrait[iInter])
          }
        }
        interactionsWithGenoTrait <- setdiff(interactionsWithGenoTrait,interactionsWithGenoToRemove)
      }else{
        interactionsWithGenoTrait <- interactionsWithGeno
      }
      ##############
      # do analysis
      if(!is.na(var(mydataSub[,"predictedValue"],na.rm=TRUE))){ # if there's variance
        if( var(mydataSub[,"predictedValue"], na.rm = TRUE) > 0 ){
          # make sure the terms to be fitted have more than one level
          if(deregress){
            mydataSub$predictedValue <- mydataSub$predictedValue/mydataSub$reliability
          }
          
          if(!is.null(randomTerm)){
            if(modelType == "rrblup"){
              reduced <- with(mydataSub,cgiarPIPE::redmm(x=designation,M=Mtrait, nPC=nPC, returnLam = TRUE))
              LGrp[["QTL"]] <- c((ncol(mydataSub)+1):(ncol(mydataSub)+ncol(reduced$Z)))
              mydataSub <- cbind(mydataSub,reduced$Z)
              rTermsTrait <- randomTerm[which(apply(data.frame(randomTerm),1,function(x){length(table(mydataSub[,x]))}) > 1)]
              rTermsTrait <- setdiff(rTermsTrait,"designation")
              rTermsTrait <- c("grp(QTL)",rTermsTrait)
            }else{
              rTermsTrait <- randomTerm[which(apply(data.frame(randomTerm),1,function(x){length(table(mydataSub[,x]))}) > 1)]
            }
          }else{rTermsTrait=NULL}
          
          if(!is.null(fixedTerm)){
            fixedTermTraitMinus <- setdiff(fixedTerm,"1") # remove intercept for a minute
            if(length(fixedTermTraitMinus) > 0){ # check that fixed effect terms have more than one level and can be fitted
              fixedTermTrait <- fixedTermTraitMinus[which(apply(data.frame(fixedTermTraitMinus),1,function(x){length(table(mydataSub[,x]))}) > 1)]
              fixedTermTrait <- c("1",fixedTermTrait)  # add it back
            }else{
              fixedTermTrait <- fixedTerm
            }
          }else{fixedTermTrait=NULL}
          
          if(!is.null(residualBy)){
            residualByTrait <- residualBy[which(apply(data.frame(residualBy),1,function(x){length(table(mydataSub[,x]))}) > 1)]
            if(length(residualByTrait) == 0){residualByTrait=NULL}
          }else{residualByTrait=NULL}
          #####
          if(length(interactionsWithGenoTrait) > 0){ # If there's interactions to be fitted build the formula terms
            if(modelType == "rrblup"){
              for(iInteraction in unique(interactionsWithGenoTrait)){ # iInteraction <- unique(interactionsWithGenoTrait)[1]
                LGrp[[paste0("QTL",iInteraction)]] <- c((ncol(mydataSub)+1):(ncol(mydataSub)+ncol(reduced$Z)))
                mydataSub <- cbind(mydataSub,reduced$Z*mydataSub[,iInteraction])
                rTermsTrait <- c(rTermsTrait,paste0("grp(QTL",iInteraction,")"))
              }
            }else{ 
              interacs <- expand.grid("designation",interactionsWithGenoTrait)
              interacs<- as.data.frame(interacs[which(as.character(interacs[,1]) != as.character(interacs[,2])),])
              interacsUnlist <- apply(interacs,1,function(x){paste(x,collapse = ":")})
              rTermsTrait <- c(rTermsTrait,interacsUnlist)
            }
          }else{
            rTermsTrait <- rTermsTrait
          }
          rTermsTrait <- setdiff(rTermsTrait, fixedTermTrait)
          
          if(length(rTermsTrait) == 0){ # if there's not a single random effect to be fitted
            ranran <- NULL
            myGinverse <- NULL # if no random effects we don't have a g inverse for LMMsolver
          }else{
            ranran <- paste("~",paste(rTermsTrait, collapse=" + ") )
          }
          fix <- paste("predictedValue ~",paste(fixedTermTrait, collapse=" + "))
          # final residual formula
          if(!is.null(residualByTrait)){ # if the user wants to fit residuals at the level of environment
            ranres <- as.formula(paste0("~",residualByTrait,""))
          }else{
            ranres <- NULL
          }
          
          mydataSub=mydataSub[with(mydataSub, order(environment)), ] # sort by environments
          mydataSub$w <- 1/(mydataSub$stdError^2) # add weights column
          if(verbose){
            cat(fix,"\n")
            cat(ranran,"\n")
          }
          
          if( modelType == "blup" ){ # if user doesn't have any special model
            # make sure the matrix only uses the leves for individuals with data
            designationFlevels <- unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"designation"])
            Ainv <- A <- diag(length(designationFlevels))
            colnames(Ainv) <- rownames(Ainv) <- designationFlevels
            colnames(A) <- rownames(A) <- designationFlevels
            inter <- character()
            onlyInA <- character() # genotypes only present in A and not in dataset
            differ <- character()
            myGinverse <- list(designation=Ainv)
            levelsInAinv <- colnames(Ainv)
            genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=designationFlevels, withMarkOnly=onlyInA)
          }else{
            if(modelType == "rrblup"){ # if user will use the rrblup 
              designationFlevels <- unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"designation"])
              inter <- character()
              onlyInA <- character() # genotypes only present in A and not in dataset
              differ <- character()
              myGinverse <- NULL
              groupTrait <- LGrp
              genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=designationFlevels, withMarkOnly=onlyInA)
            }else{ # if user wants to do a gblup, pblup or ssgblup model
              designationFlevels <- as.character(unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"designation"]))
              
              if(modelType %in% c("pblup","ssgblup")){ # we need to calculate NRM
                N <- cgiarBase::nrm2(pedData=phenoDTfile$data$pedigree)
              }
              if(modelType %in% c("gblup","ssgblup")){ # we need to calculate GRM
                commonBetweenMandP <- intersect(rownames(phenoDTfile$data$geno),designationFlevels)
                if(length(commonBetweenMandP) < 2){ stop("Markers could not be matched with phenotypes. Please ensure that you have used the right marker file or check the rownames of your marker matrix and ids of your phenotypes.", call. = FALSE)                }
                M <- phenoDTfile$data$geno[commonBetweenMandP,]
                if(ncol(M) > 5000){ # we remove that many markers if a big snp chip
                  A <- sommer::A.mat(M[,sample(1:ncol(M), 5000)])
                }else{ A <- sommer::A.mat(M) };  M <- NULL
                if(modelType == "ssgblup"){ # only if ssgblup we merge
                  A <- sommer::H.mat(N,A, tau=1,  omega=1, tolparinv=1e-6)
                }
              }
              if(modelType %in% c("pblup")){ A <- N  }
              
              badGeno <- which(rownames(A) == "") # should have no covariance with anyone
              if(length(badGeno) > 0){A[badGeno,2:ncol(A)]=0; A[2:nrow(A),badGeno]=0} # make zero covariance with this genotype
              badBlankGenotype <- which(colnames(A)=="")
              if(length(badBlankGenotype) > 0){A <- A[-badBlankGenotype,-badBlankGenotype]}
              inter <- intersect(designationFlevels,colnames(A)) # go for sure
              onlyInA <- setdiff(rownames(phenoDTfile$data$geno),designationFlevels)
              differ <- setdiff(designationFlevels,inter) # are missing in A matrix
              genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=differ, withMarkOnly=onlyInA)
              # get inverse matrix
              if(length(inter) > 0){ #
                A1inv <- solve(A[inter,inter] + diag(tolParInv,length(inter), length(inter)))
              }else{A1inv <- matrix(0,0,0)}
              if(length(differ) > 0){ # we have to add individuals without markers or not being part of the GRM?
                if(length(differ) > 1){ # there's at least 2 inds to be added
                  A2inv <- diag(x=rep(mean(diag(A1inv)),length(differ)))
                }else{ A2inv <- diag(1)*mean(diag(A1inv)) } # there's only one individual to be added
                colnames(A2inv) <- rownames(A2inv) <- differ
              }else{A2inv <- matrix(0,0,0)}
              Ainv <- sommer::adiag1(A1inv,A2inv)
              Ainv[lower.tri(Ainv)] <- t(Ainv)[lower.tri(Ainv)] # fill the lower triangular
              colnames(Ainv) <- rownames(Ainv) <- c(colnames(A1inv), colnames(A2inv))
              A1inv <- NULL; A2inv <- NULL;
              levelsInAinv <- colnames(Ainv)
              myGinverse <- list(designation=Ainv)
            }
          }
          if(length(ranran) == 0){ranFormulation=NULL}else{ranFormulation=as.formula(ranran)}
          if(useWeights){
            weightsFormulation="w"
            if(verbose){print("Using weights in the analysis. Residual variance will be fixed to 1.")  }
          }else{
            weightsFormulation=NULL
            if(verbose){print("Ignoring weights in the analysis. Residual variance will be estimated.")  }
          }
          mix <- try(
            LMMsolver::LMMsolve(fixed =as.formula(fix),
                                random = ranFormulation,
                                residual=ranres,
                                weights = weightsFormulation,
                                ginverse = myGinverse,
                                group = groupTrait,
                                family = eval(parse(text = traitFamily[iTrait])),
                                data = mydataSub, maxit = maxIters),
            silent = TRUE
          )
          myGinverse <- NULL      #
          currentModeling <- data.frame(module="mta", analysisId=mtaAnalysisId,trait=iTrait, environment="across", parameter=c("fixedFormula","randomFormula","family"), value=c(fix,ranran,traitFamily[iTrait]))
          phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
          # print(mix$VarDf)
          if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model
            if(is.null(phenoDTfile$metadata$weather)){numericMetas <- character()}
            for(iIndex in c(numericMetas,"envIndex")){
              if( (iIndex %in% interactionsWithGenoTrait) ){names(mix$ndxCoefficients[[paste0("designation:",iIndex)]]) <- names(mix$ndxCoefficients$designation) } # copy same names than main designation effect
            }
            
            iGenoUnit <- "designation" # in MET iGenoUnit is always "designation" only in STA we allow for different
            
            ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(modelType == "rrblup"){ # user wants to predict using marker effects
              ##
              Mfull <- phenoDTfile$data$geno
              if((ncol(Mfull) < nrow(Mfull)) | nPC==0){M2 <- Mfull}else{  M2 <- tcrossprod(Mfull)}
              xx2 = with(mydataSub, cgiarPIPE::redmm(x=designation, M=M2, nPC=nPC, returnLam = TRUE)) # we need the new lambda for the fullt marker matrix
              ss <- mix$VarDf;  rownames(ss) <- ss$VarComp
              pp <- list()
              for(iGroup in names(LGrp)){ # iGroup <- names(LGrp)[1]  for each rrBLUP effect
                shouldBeOne <- which(mix$ndxCoefficients[[iGroup]] == 0)
                if(length(shouldBeOne) > 0){mix$ndxCoefficients[[iGroup]][shouldBeOne] = 1}
                blup <- mix$coefMME[mix$ndxCoefficients[[iGroup]]]
                names(blup) <- names(mix$ndxCoefficients[[iGroup]])
                names(blup) <- gsub(paste0(iGroup,"_"),"",names(blup))
                predictedValue <- (xx2$Lam %*% blup[colnames(xx2$Lam)]) + mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
                if(length(shouldBeOne) > 0){predictedValue[1] = mix$coefMME[mix$ndxCoefficients$`(Intercept)`]}
                # names(predictedValue) <- names(mix$ndxCoefficients[[iGroup]])
                dims <- mix$EDdf
                start <- sum(dims[1:(which(dims$Term == iGroup) - 1),"Model"]) # we don't add a one because we need the intercept
                nEffects <- length(mix$coefMME[mix$ndxCoefficients[[iGroup]]])
                pev <- as.matrix(solve(mix$C))[start:(start+nEffects-1),start:(start+nEffects-1)]
                pev <- xx2$Lam %*% pev %*% t(xx2$Lam) 
                stdError <- (sqrt(Matrix::diag(pev)))
                reliability <- 1 - (stdError/as.numeric(var(predictedValue)))
                badRels <- which(reliability > 1); if(length(badRels) > 0){reliability[badRels] <- 0.9999}
                badRels2 <- which(reliability < 0); if(length(badRels2) > 0){reliability[badRels2] <- 0}
                pp[[iGroup]] <- data.frame(designation=rownames(predictedValue), predictedValue=predictedValue, stdError=stdError, reliability=reliability,
                                           trait=paste0(iTrait,"_",iGroup) )
                ## save metrics
                phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                             data.frame(module="mta",analysisId=mtaAnalysisId, trait= paste0(iTrait,"_",iGroup), environment="across", 
                                                        parameter=c("mean","CV", "r2","Vg"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML"), 
                                                        value=c(mean(pp$predictedValue, na.rm=TRUE), cv, median(pp$reliability), ss[iGroup,"Variance"]), 
                                                        stdError=c(0,0,sd(pp$reliability, na.rm = TRUE)/sqrt(length(pp$reliability)),0)
                                             )
                )
                counter <- counter+1
              }
              pp <- do.call(rbind,pp)
              ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            }else{ # user wants to do genetic evaluation
              
              
              if(length(grep("designation", ranran)) > 0){ # geno was random
                ######################################################
                ## do genetic evaluation only if there was genotype as random effect
                ######################################################
                if(modelType %in% c("gblup","ssgblup")){ # if user provided markers
                  allIndsNames <- rownames(phenoDTfile$data$geno) # individuals in the marker matrix
                }else{# if user provided relationship matrix
                  allIndsNames <- colnames(A) # individuals in the relationship matrix
                }
                indNamesInAinv <- levelsInAinv # reduced set of individuals in the met analysis
                indNamesInAinvNoDiffer <- setdiff(indNamesInAinv,differ) # individuals with no markers data that were just expanded in the A matrix
                indsToBePredicted <- c(indNamesInAinv, setdiff(allIndsNames, indNamesInAinv) )# additional individuals to be predicted in the A matrix
                # Ainv <- NULL
                if(length(indsToBePredicted) > 0 ){ # if there's additional inds to be predicted
                  groups0 <- sort(rep(paste0("g",1:1000),batchSizeToPredict)) # groups of 500
                  grouping <- data.frame(id=indsToBePredicted,grp=groups0[1:length(indsToBePredicted)])
                  grouping$grp[1:length(indNamesInAinv)] <- "g1" ## all the ones in the model go to group 1
                  groups <- unique(grouping$grp)
                }else{ # no additional inds
                  grouping <- data.frame(id=indNamesInAinv,grp="g1")
                  groups <- unique(grouping$grp)
                }
                nEffectsToStore <- length(interactionsWithGenoTrait)+1
                genEva <- lapply(vector(mode="list",nEffectsToStore), function(x,nn){vector(mode="list",nn)}, nn=4) #list of lists
                names(genEva) <- names(mix$ndxCoefficients)[grep("designation",names(mix$ndxCoefficients))]
                genEva <- lapply(genEva, function(x){names(x) <- c("predictedValue","stdError","pev","R2"); return(x)})
                ttGrouping <- table(grouping$grp)
                for(iGroup in 1:length(groups)){ # iGroup=1  # for each group do the genetic evaluation
                  if(verbose){print(paste("Predicting batch",iGroup))}
                  iGrupNames <- grouping[which(grouping$grp == groups[iGroup]),"id"]
                  # make a new A for the ith group
                  present <- unique(c( intersect(colnames(A), indNamesInAinv), iGrupNames ))
                  mydataSub2 <- mydataSub[which(!is.na(mydataSub$predictedValue)),]
                  mydataSub2 <- droplevels(mydataSub2[which(mydataSub2$designation %in% c(present)),])
                  if(modelType %in% c("gblup","ssgblup")){ # if user provided markers
                    AGE <- sommer::A.mat(phenoDTfile$data$geno[intersect(present,rownames(phenoDTfile$data$geno)),])
                  }else{
                    present2 <- intersect(colnames(A),present)
                    AGE <- A[present2,present2]
                  }
                  ######################################
                  ## do the genetic evaluation
                  genEvaProv <- cgiarBase::geneticEvaluation3(
                    fixed =as.formula(fix),
                    random = as.formula(ranran),
                    rcov=ranres,
                    weights = "w",
                    data = mydataSub2,
                    varComp=mix$VarDf,
                    keep=ttGrouping[iGroup],
                    AFGE=AGE
                  )
                  AGE <- NULL
                  # print("finished")
                  # fill the list
                  for(uu in 1:length(genEvaProv)){ # for each genetic effect
                    genEva[[uu]]$predictedValue <- c(genEva[[uu]]$predictedValue,genEvaProv[[uu]]$predictedValue )
                    genEva[[uu]]$stdError <- c(genEva[[uu]]$stdError,genEvaProv[[uu]]$stdError )
                    if(iGroup == 1){
                      genEva[[uu]]$pev <- genEvaProv[[uu]]$pev
                      genEva[[uu]]$R2 <- genEvaProv[[uu]]$R2
                    }else{
                      genEva[[uu]]$pev <- Matrix::bdiag(genEva[[uu]]$pev,genEvaProv[[uu]]$pev )
                      genEva[[uu]]$R2 <- c(genEva[[uu]]$R2,genEvaProv[[uu]]$R2 ) #Matrix::bdiag(genEva[[uu]]$R2,genEvaProv[[uu]]$R2 )
                    }
                  }
                  genEvaProv <- NULL
                }
                ######################################################
                ## end of genetic evaluation function
                ######################################################
              }
              ## now extract the needed values
              if(iGenoUnit %in% fixedTermTrait){ # user wants fixed effect predictions for genotype
                shouldBeOne <- which(mix$ndxCoefficients[[iGenoUnit]] == 0)
                if(length(shouldBeOne) > 0){mix$ndxCoefficients[[iGenoUnit]][shouldBeOne] = 1}
                predictedValue <- mix$coefMME[mix$ndxCoefficients[[iGenoUnit]]] + mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
                if(length(shouldBeOne) > 0){predictedValue[1] = mix$coefMME[mix$ndxCoefficients$`(Intercept)`]}
                names(predictedValue) <- names(mix$ndxCoefficients[[iGenoUnit]])
                dims <- mix$EDdf
                start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) # we don't add a one because we need the intercept
                pev <- as.matrix(solve(mix$C))[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                stdError <- (sqrt(diag(pev)))
              }else{ # user wants random effect predictions for genotype (main effect)
                predictedValue <- genEva$designation$predictedValue
                stdError <- genEva$designation$stdError
                pev <- genEva$designation$pev
              }
              designation <- gsub("designation","", names(predictedValue))
              pp <- data.frame(designation,predictedValue,stdError)
              ss = mix$VarDf; rownames(ss) <- ss$VarComp
              Vg <- ss["designation",2]; Vr <- ss["residual",2]
              if(iGenoUnit %in% fixedTermTrait){ # add reliabilities to the data frame
                pp$reliability <- NA
              }else{ # if random, reliability can be calculated for main effect
                # print(names(genEva))
                pp$reliability <- genEva$designation$R2
                badRels <- which(pp$reliability > 1); if(length(badRels) > 0){pp$reliability[badRels] <- 0.9999}
                badRels2 <- which(pp$reliability < 0); if(length(badRels2) > 0){pp$reliability[badRels2] <- 0}
              }
              pp$trait <- iTrait # add trait
              cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
              ## save metrics
              phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                           data.frame(module="mta",analysisId=mtaAnalysisId, trait=iTrait, environment="across", 
                                                      parameter=c("mean","CV", "r2","Vg"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML"), 
                                                      value=c(mean(pp$predictedValue, na.rm=TRUE), cv, mean(pp$reliability), Vg), 
                                                      stdError=c(0,0,sd(pp$reliability, na.rm = TRUE)/sqrt(length(pp$reliability)),0)
                                           )
              )
              ## genetic variances
              lpv <- sum(mix$EDdf$Model[1:which(mix$EDdf$Term == "designation")])+1 # to be used as a starting point if random regression is requested
              # extract sensitivities if interaction is requested
              if(length(interactionsWithGenoTrait) > 0){ # if there's interactions
                if( length(intersect(interactionsWithGenoTrait, c("envIndex","timePoint","latitude","longitude","altitude","weather1","weather2"))) > 0 ){
                  
                  for(iInteractionTrait in c("envIndex","timePoint","latitude","longitude","altitude","weather1","weather2")){ # iInteractionTrait="envIndex"
                    
                    if( (iInteractionTrait %in% interactionsWithGenoTrait) ){ # iInteractionTrait="envIndex"
                      
                      counter <- counter+1
                      iGenoUnit <- paste0("designation:",iInteractionTrait)#"designation:envIndex"
                      if(iGenoUnit %in% fixedTermTrait){ # user wants fixed effect predictions for genotype:envIndex
                        shouldBeOne <- which(mix$ndxCoefficients[[iGenoUnit]] == 0)
                        if(length(shouldBeOne) > 0){mix$ndxCoefficients[[iGenoUnit]][shouldBeOne] = 1}
                        predictedValue <- mix$coefMME[mix$ndxCoefficients[[iGenoUnit]]] #+ mix$coefficients$`(Intercept)`
                        dims <- mix$EDdf
                        start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) # we don't add a one because we need the intercept
                        stdError <- (sqrt(diag(as.matrix(solve(mix$C)))))[start:(start+length(predictedValue)-1)]
                      }else{ # user wants random effect predictions for genotype:envIndex
                        predictedValue <- genEva[[iGenoUnit]]$predictedValue #+ mix$coefficients$`(Intercept)`
                        stdError <- genEva[[iGenoUnit]]$stdError
                      }
                      designation <- gsub("designation","", names(predictedValue))
                      pp2 <- data.frame(designation,predictedValue,stdError)
                      Vg <- ss[iGenoUnit,2];
                      if(iGenoUnit %in% fixedTermTrait){pp$reliability <- 1e-6}else{pp2$reliability <- genEva[[iGenoUnit]]$R2}
                      pp2$trait <- paste(iTrait,iInteractionTrait,sep="-")
                      cv <- (sd(pp2$predictedValue,na.rm=TRUE)/mean(pp2$predictedValue,na.rm=TRUE))*100
                      ## save metrics
                      phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                                   data.frame(module="mta",analysisId=mtaAnalysisId, trait=paste(iTrait,iInteractionTrait,sep="-"), 
                                                              environment="across", 
                                                              parameter=c("mean","CV", "r2","Vg"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML"), 
                                                              value=c(mean(pp2$predictedValue, na.rm=TRUE), cv, mean(pp2$reliability), Vg), 
                                                              stdError=c(0,0,sd(pp2$reliability, na.rm = TRUE)/sqrt(length(pp2$reliability)),0)
                                                   )
                      )
                      pp <- rbind(pp,pp2) # bind predictions
                      
                    }
                    
                  }
                  
                }
              } # if there's even interactions
              
            }
            ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
          }else{ # if model failed
            if(modelType == "rrblup"){
              # do nothing
            }else{
              if(verbose){ cat(paste("Mixed model failed for this combination. Aggregating and assuming h2 = 0 \n"))}
              pp <- aggregate(predictedValue ~ designation, FUN=mean, data=mydataSub)
              pp$stdError <- sd(pp$predictedValue)  # 1
              # pp$status <- "Aggregated"
              pp$reliability <- 1e-6
              pp$trait <- iTrait
              cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
              ## save metrics
              phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                           data.frame(module="mta",analysisId=mtaAnalysisId, trait=iTrait, 
                                                      environment="across", 
                                                      parameter=c("mean","CV", "r2","Vg"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML"), 
                                                      value=c(mean(pp$predictedValue, na.rm=TRUE), cv, 0, 0), 
                                                      stdError=c(0,0,0,0)
                                           )
              )
            }
            
          }
          
          pp$environment <- "across"
          mydataForEntryType <- droplevels(mydata[which(mydata$trait == iTrait),])
          pp$entryType <- apply(data.frame(pp$designation),1,function(x){
            found <- which(mydataForEntryType$designation %in% x)
            if(length(found) > 0){
              x2 <- paste(sort(unique(toupper(trimws(mydataForEntryType[found,"entryType"])))), collapse = "#");
            }else{x2 <- "unknown"}
            return(x2)
          })
          mydataForEntryType <- NULL
          if(modelType == "rrblup"){
            pp$entryType <- "markerEffect"
          }else{
            pp$entryType <- paste(ifelse(as.character(pp$designation) %in% differ, "TGV", surrogate),
                                  pp$entryType,
                                  ifelse(as.character(pp$designation) %in% onlyInA, "predicted", "tested"),
                                  sep="_")
          }
          ###
          predictionsList[[counter2]] <- pp;
          counter=counter+1
        }
      }
    }
    counter2 = counter2+1
  }
  #
  trait <- setdiff(trait,traitToRemove)
  heritLB <- heritLB[which(traitOrig %in% trait)]
  heritUB <- heritUB[which(traitOrig %in% trait)]
  if(length(predictionsList) == 0){stop("There was no predictions to work with. Please look at your H2 boundaries. You may be discarding all envs.",call. = FALSE)}
  predictionsBind <- do.call(rbind, predictionsList)
  predictionsBind$analysisId <- mtaAnalysisId
  ##########################################
  ## add timePoint of origin, stage and designation code
  entries <- unique(mydata[,"designation"])
  baseOrigin <- do.call(rbind, apply(data.frame(entries),1,function(x){
    out1 <- (sort(mydata[which(mydata$designation %in% x),"gid"], decreasing = FALSE))[1]
    out2 <- (sort(mydata[which(mydata$designation %in% x),"mother"], decreasing = FALSE))[1]
    out3 <- (sort(mydata[which(mydata$designation %in% x),"father"], decreasing = FALSE))[1]
    out4 <- paste(unique(sort(mydata[which(mydata$designation %in% x),"pipeline"], decreasing = FALSE)),collapse=", ")
    y <- data.frame(designation=x,gid=out1,mother=out2,father=out3,pipeline=out4)
    return(y)
  }))
  predictionsBind <- merge(predictionsBind,baseOrigin, by="designation", all.x=TRUE)
  predictionsBind$module <- "mta"
  #########################################
  ## update databases
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions,
                                   predictionsBind[,colnames(phenoDTfile$predictions)])
  phenoDTfile$status <- rbind( phenoDTfile$status, data.frame(module="mta", analysisId=mtaAnalysisId))
  ## add which data was used as input
  modeling <- data.frame(module="mta",  analysisId=mtaAnalysisId, trait=c("inputObject","all"), environment="general", 
                         parameter= c("analysisId","estimateType"), value= c(analysisId,ifelse("designation"%in% randomTerm,"random","fixed") ))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)
}
