###Generate the potential Z Design Matrix for the four potential models

#' Title
#'
#' @param tumor.character.cat
#' @param tumor.number
#' @param tumor.names
#' @param freq.subtypes
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param tumor.character.cat
#' @param tumor.number
#' @param tumor.names
#' @param freq.subtypes
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param tumor.names
#'
#' @return
#' @export
#'
#' @examples
GenerateZDesignNamesMainEffect <- function(tumor.names){

  z.design.names.main.effect <- "baseline effect"

  z.design.names.main.effect <- c(z.design.names.main.effect,
                                  paste0(tumor.names," main effect"))

  return(z.design.names.main.effect)

}

###Generate z design matrix for pairwise interaction model
#' Title
#'
#' @param tumor.character.cat
#' @param tumor.number
#' @param tumor.names
#' @param freq.subtypes
#'
#' @return
#' @export
#'
#' @examples
GenerateZDesignPairwiseInteraction <- function(tumor.character.cat,tumor.number,tumor.names,freq.subtypes){
  cutoff <- 10
  z.design.pairwise.interaction <-
    GenerateZDesignMainEffect(tumor.character.cat,
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
#' Title
#'
#' @param tumor.character.cat
#' @param tumor.number
#' @param tumor.names
#' @param freq.subtypes
#'
#' @return
#' @export
#'
#' @examples
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
