##z design transform to z all
ZDesigntoZall <- function(baselineonly,
                          main.ffect,
                          pairwise.interaction,
                          saturated,
                          z.design.baselineonly,
                          z.design.main.effect,
                          z.design.pairwise.interaction,
                          z.design.saturated) {
  M <- nrow(z.design.main.effect)
  ##the number of covariates in different potential model structures
  baselineonly.number <- CountCovarNumber(baselineonly)
  main.effect.number <- CountCovarNumber(main.effect)
  pairwise.interaction.number <- CountCovarNumber(pairwise.interaction)
  saturated.number <- CountCovarNumber(saturated)
  ###second.stage.category for different model structures
  baselineonly.second.cat <- 1
  main.effect.second.cat <- ncol(z.design.main.effect)
  pairwise.interaction.second.cat <- ncol(z.design.pairwise.interaction)
  saturated.second.cat <- ncol(z.design.saturated)
  ###1 for intercept
  total.covar.number <- 1+ baselineonly.number+main.effect.number+
    pairwise.interaction.number+saturated.number

  z.all <- matrix(0,nrow=(M*total.covar.number),ncol = (M+baselineonly.number*baselineonly.second.cat+
                                                          main.effect.second.cat*main.effect.number+
                                                          pairwise.interaction.second.cat*pairwise.interaction.number)+saturated.second.cat*saturated.number)

  for(i in c("intercept","baselineonly",
             "main.effect","pairwise.interection",
             "satuared")){
    ##we always keep intercept as saturated model and to simply, we always use diagnonal matrix for intercept
    if(i=="intercept"){
      ###row start and column start point for this category
      row.start <- 0
      column.start <- 0
      for(j in 1:M){
        z.all[row.start+1+(j-1)*total.covar.number,(column.start+j)] = 1
      }
    }else if(i=="baselineonly"){
      column.start = M
      row.start <- 1
      ###test whether there is any baselineonly variable
      if(baselineonly.number!=0){
        for(j in 1:M){
          for(k in 1:baseline.number){
            z.all[row.start+k+(j-1)*total.covar.number,(column.start+k)] <- 1
          }
        }
      }
    }else if(i=="main.effect"){
      column.start <- M+baselineonly.number
      row.start <- 1+baselineonly.number
      if(main.effect.number!=0){
        for(j in 1:M){
          for(k in 1:main.effect.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*main.effect.second.cat+1):
                    (column.start+k*main.effect.second.cat)] <- as.matrix(z.design.main.effect[j,])
          }
        }
      }
    }else if(i == "pairwise.interaction"){
      column.start <- M+baselineonly.number+main.effect.number*main.effect.second.cat
      row.start <- 1+baselineonly.number+main.effect.number
      if(pairwise.interaction.number!=0){
        for(j in 1:M){
          for(k in 1:pairwise.interaction.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*pairwise.interaction.second.cat+1):
                    (column.start+k*pairwise.interaction.second.cat)] <- as.matrix(z.design.pairwise.interaction[j,])
          }
        }
      }
    }else {
      column.start <- M+baselineonly.number+main.effect.number*main.effect.second.cat+
        pairwise.interaction.number*pairwise.interaction.second.cat
      row.start <- 1+baselineonly.number+main.effect.number+pairwise.interaction.number
      if(saturated.number!=0){
        for(j in 1:M){
          for(k in 1:saturated.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*saturated.second.cat+1):
                    (column.start+k*saturated.second.cat)] <- as.matrix(z.design.saturated[j,])
          }
        }
      }
    }



  }
  return(z.all)
}

##count variable number
CountCovarNumber <- function(covar) {
  if(is.null(covar)){
    return(0)
  }else if(is.vector(covar)){
    return(1)
  }else{
    return(ncol(covar))
  }
}

##check whether control subtypes tumor characteristics is NA
CheckControlTumor <- function(y.case.control,y.tumor){
  idx <- which(y.case.control==0)
  y.tumor.control = y.tumor[idx,]
  return(any(!is.na(y.tumor.control)))


}



###Generate the complete y pheno dataframe based on the incomplete data file
GenerateCompleteYPheno <- function(y.pheno,missingTumorIndicator){
  tumor.number  <-  ncol(y.pheno)-1
  find.missing.position.text = "idx <- which("
  for(i in 2:(tumor.number+1)){
    if(i == (tumor.number+1)){
      find.missing.position.text <- paste0(find.missing.position.text,"y.pheno[,",i,"]==missingTumorIndicator)")
    }else{
      find.missing.position.text <- paste0(find.missing.position.text,"y.pheno[,",i,"]==missingTumorIndicator|")
    }
  }
  eval(parse(text=find.missing.position.text))
  if(length(idx)!=0){
    y.pheno.complete = y.pheno[-idx,]
  }else{
    y.pheno.complete = y.pheno
  }
  return(y.pheno.complete)
}

### Generate the frequency of subtypes in the complete data file
GenerateFreqTable <- function(y.pheno.complete){
  eval_text = NULL
  tumor.number = ncol(y.pheno.complete)-1
  for(i in 2:(tumor.number+1)){
    if(i==(tumor.number+1)){
      eval_text = paste0(eval_text,paste0("y.pheno.complete[,",i,"]"))
    }else{
      eval_text = paste0(eval_text,paste0("y.pheno.complete[,",i,"],"))
    }
  }
  eval_text = paste0("as.data.frame(table(",eval_text,"),stringsAsFactors = F)")
  result = eval(parse(text=eval_text))
  result <- apply(result,2,as.numeric)
  return(result)
}
### Generate the potential tumor characteristic category(binary or categorical)
GenerateTumorCharacterCat <- function(y.pheno.complete){
  tumor.number = ncol(y.pheno.complete)-1
  y.tumor.complete = y.pheno.complete[,2:(tumor.number+1)]
  tumor.character.cat = list()
  for(i in 1:tumor.number){
    unique.tumor.cat = unique(y.tumor.complete[!is.na(y.tumor.complete[,i]),i])
    unique.tumor.cat = unique.tumor.cat[order(unique.tumor.cat)]
    tumor.character.cat[[i]] = unique.tumor.cat
  }
  return(tumor.character.cat)
}

###Generate the potential Z Design Matrix for the four potential models
GenerateZDesignBaselineonly <- function(tumor.character.cat,tumor.number,tumor.names,freq.subtypes){
  M = 1
  cutoff <- 10
  for(i in 1:tumor.number){
    M = M*length(tumor.character.cat[[i]])
  }

  z.design.baselineonly <- matrix(rep(1,M),M,1)

  freq = freq.subtypes[,ncol(freq.subtypes)]
  idx <- which(freq<=cutoff)
  if(length(idx!=0)){
    return(z.design.baselineonly[-idx,])
  }else{
    return(z.design.baselineonly)
  }

}

###Generate z design matrix for main effect model
GenerateZDesignMainEffect <- function(tumor.character.cat,tumor.number,tumor.names,freq.subtypes){
  z.design.main.effect.text <- NULL
  cutoff <- 10
  for(i in 1:tumor.number){
    if(i==tumor.number){
      z.design.main.effect.text <- paste0(z.design.main.effect.text,
                                          "tumor.character.cat[[",i,"]]")
    }else{
      z.design.main.effect.text <- paste0(z.design.main.effect.text,
                                          "tumor.character.cat[[",i,"]],")
    }
  }
  z.design.main.effect.text <- paste0("z.design.main.effect <- expand.grid(",
                                      z.design.main.effect.text,
                                      ")")
  eval(parse(text=z.design.main.effect.text))
  z.design.main.effect <- cbind(1,z.design.main.effect)
  colnames(z.design.main.effect) <- GenerateZDesignNamesMainEffect(tumor.names)
  freq = freq.subtypes[,ncol(freq.subtypes)]
  idx <- which(freq<=cutoff)
  if(length(idx!=0)){
    return(z.design.main.effect[-idx,])
  }else{
    return(z.design.main.effect)
  }
}

###Generate tumor.names for z design matrix for main effect
GenerateZDesignNamesMainEffect <- function(tumor.names){

  z.design.names.main.effect <- "baseline effect"

  z.design.names.main.effect <- c(z.design.names.main.effect,
                                  paste0(tumor.names," main effect"))

  return(z.design.names.main.effect)

}

###Generate z design matrix for pairwise interaction model
GenerateZDesignPairwiseInteraction <- function(tumor.character.cat,tumor.number,tumor.names,freq.subtypes){
  cutoff <- 10
  z.design.pairwise.interaction <- GenerateZDesignMainEffect(tumor.character.cat,
                                                             tumor.number,
                                                             tumor.names,
                                                             freq.subtypes)
  z.design.names.pairwise.interaction <- colnames(z.design.pairwise.interaction)
  all.pairwise.combnation <- combn(tumor.number,2)+1
  combn.number <- ncol(all.pairwise.combnation)

  for(i in 1:combn.number){
    col1 <- all.pairwise.combnation[1,i]
    col2 <- all.pairwise.combnation[2,i]
    newcol <-  z.design.pairwise.interaction[,col1]*
      z.design.pairwise.interaction[,col2]
    z.design.names.pairwise.interaction <- c(z.design.names.pairwise.interaction,
                                             #col1-1 is due to there is basline effect in the z design matrix
                                             paste0(tumor.names[col1-1],
                                                    tumor.names[col2-1],
                                                    " interaction effect"))
    z.design.pairwise.interaction <- cbind(z.design.pairwise.interaction,
                                           newcol)
  }
  colnames(z.design.pairwise.interaction) <- z.design.names.pairwise.interaction
  freq = freq.subtypes[,ncol(freq.subtypes)]
  idx <- which(freq<=cutoff)
  if(length(idx!=0)){
    return(z.design.pairwise.interaction[-idx,])
  }else{
    return(z.design.pairwise.interaction)
  }



}
##Generate the z design matrix for saturated model
GenerateZDesignSaturated <- function(tumor.character.cat,tumor.number,tumor.names,freq.subtypes) {
  cutoff <- 10
  z.design.saturated <- GenerateZDesignMainEffect(tumor.character.cat,
                                                  tumor.number,
                                                  tumor.names,
                                                  freq.subtypes)
  z.design.names.saturated <- colnames(z.design.saturated)
  ##j represent the order of the interaction
  for(j in 2:tumor.number){
    all.combnation <- combn(tumor.number,j)+1
    combn.number <- ncol(all.combnation)
    #####i represent the ith combination within the jth order interaction
    for(i in 1:combn.number){
      newcol <- rep(1,nrow(z.design.saturated))
      ###k present the kth tumor characteristics within the ith combination
      z.design.names.saturated.one.column <- NULL
      for(k in 1:j){
        col.number <- all.combnation[k,i]
        newcol <-  newcol*z.design.saturated[,col.number]
        z.design.names.saturated.one.column <- paste0(z.design.names.saturated.one.column
                                                      ,tumor.names[col.number-1])
      }
      z.design.names.saturated.one.column <- paste0(z.design.names.saturated.one.column,
                                                    " interaction effect")
      z.design.names.saturated <- c(z.design.names.saturated,
                                    z.design.names.saturated.one.column)
      z.design.saturated <- cbind(z.design.saturated,newcol)

    }

  }
  colnames(z.design.saturated) <- z.design.names.saturated

  freq = freq.subtypes[,ncol(freq.subtypes)]
  idx <- which(freq<=cutoff)
  if(length(idx!=0)){
    return(z.design.saturated[-idx,])
  }else{
    return(z.design.saturated)
  }

}


####set up the start vaue for the parameters
StartValueFunction = function(freq.subtypes,y.case.control,z.all){
  ###cutoff for take one subject
  cutoff=10
  ncontrol <- sum(y.case.control==0)
  p.freq <- ncol(freq.subtypes)
  freq = freq.subtypes[,p.freq]
  idx =which(freq<=cutoff)
  if(length(idx)!=0){
    freq.subtypes = freq.subtypes[-idx,]
    freq = freq.subtypes[,p.freq]
    total = sum(freq)+ncontrol
    p.empirical = freq/ncontrol
    delta_inter = log(p.empirical/(1-p.empirical))
    return(list(delta_inter,idx))
  }else{
    total = sum(freq)+ncontrol
    p.empirical = freq/ncontrol
    delta_inter = log(p.empirical/(1-p.empirical))
  }
  delta0 <- rep(0,ncol(z.all))
  delta0[1:length(delta_inter)] <- delta_inter
  return(delta0)
}
