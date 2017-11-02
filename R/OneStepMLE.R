
#' Title
#'
#' @param y
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param missingTumorIndicator
#'
#' @return
#' @export
#'
#' @examples
OneStepMLE <- function(y,
                     baselineonly,
                     additive,
                     pairwise.interaction,
                     saturated,
                     missingTumorIndicator){

  tumor.number <- ncol(y)-1
  y.case.control <- y[,1]
  y.tumor <- y[,2:(tumor.number+1)]
  y.pheno.complete <- GenerateCompleteYPheno(y,missingTumorIndicator)
  freq.subtypes <- GenerateFreqTable(y.pheno.complete)
  if(CheckControlTumor(y.case.control,y.tumor)==1){
    return(print("ERROR:The tumor characteristics for control subtypes should put as NA"))
  }
  tumor.names <- colnames(y.tumor)
  if(is.null(tumor.names)){
    tumor.names <- paste0(c(1:tumor.number))
  }
  tumor.character.cat = GenerateTumorCharacterCat(y.pheno.complete)







  z.design.baselineonly <- GenerateZDesignBaselineonly(tumor.character.cat,
                                                       tumor.number,
                                                       tumor.names,
                                                       freq.subtypes)
  z.design.additive <- GenerateZDesignAdditive(tumor.character.cat,
                                               tumor.number,
                                               tumor.names,
                                               freq.subtypes)
  z.design.pairwise.interaction <- GenerateZDesignPairwiseInteraction(tumor.character.cat,
                                                                      tumor.number,
                                                                      tumor.names,
                                                                      freq.subtypes)
  z.design.saturated <- GenerateZDesignSaturated(tumor.character.cat,
                                                 tumor.number,
                                                 tumor.names,
                                                 freq.subtypes)
  full.second.stage.names <- colnames(z.design.saturated)
  covar.names <- GenerateCovarName(baselineonly,
                                   additive,
                                   pairwise.interaction,
                                   saturated)

  z.all <- ZDesigntoZall(baselineonly,
                         additive,
                         pairwise.interaction,
                         saturated,
                         z.design.baselineonly,
                         z.design.additive,
                         z.design.pairwise.interaction,
                         z.design.saturated)
  delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
  #x.all has no intercept yet
  #we will add the intercept in C code
  x.all <- GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated)
  ###z standard matrix means the additive model z design matrix without baseline effect
  ###z standard matrix is used to match the missing tumor characteristics to the complete subtypes

  y <- as.matrix(y)
  x.all <- as.matrix(x.all)
  z.standard <- z.design.additive[,-1]
  M <- as.integer(nrow(z.standard))
  p.main <- ncol(z.standard)+1

  EM.result = EMStep(delta0,as.matrix(y),x.all,z.standard,z.all,missingTumorIndicator)
  ###delta represent second stage parameters
  delta <- EM.result$delta
  covariance.delta <- solve(EM.result$infor_obs)
  loglikelihood <- EM.result$loglikelihood
  AIC <- EM.result$AIC
  second.stage.mat <-
    GenerateSecondStageMat(baselineonly,
                           additive,
                           pairwise.interaction,
                           saturated,
                           M,
                           full.second.stage.names,
                           covar.names,
                           delta,
                           z.design.additive,
                           z.design.pairwise.interaction,
                           z.design.saturated)
  ##take out the intercept from second stage parameters

  takeout.intercept.result <- TakeoutIntercept(delta,covariance.delta,
                                               M,
                                               tumor.names,
                                               z.all,covar.names)
  beta <- takeout.intercept.result$beta
  covariance.beta <- takeout.intercept.result$covariance.beta
  delta.no.inter <- takeout.intercept.result$delta.no.inter
  covariance.delta.no.inter <-
    takeout.intercept.result$covariance.delta.no.inter
  beta.no.inter <- takeout.intercept.result$beta.no.inter
  covariance.beta.no.inter <- takeout.intercept.result$covariance.beta.no.inter



  second.stage.test <- SecondStageTest(delta.no.inter,covariance.delta.no.inter,M,second.stage.mat)
  global.test <- GenerateGlobalTest(delta.no.inter,
                                    covariance.delta.no.inter,
                                    M,
                                    second.stage.mat)
  ##beta represent first stage parameters

  subtypes.names <- GenerateSubtypesName(z.design.additive,M,
                                         tumor.names)
  first.stage.mat <- GenerateFirstStageMat(beta,
                                           covar.names,
                                           subtypes.names)

  first.stage.test <- FirstStageTest(beta.no.inter,
                                     covariance.beta.no.inter,
                                     M,
                                     first.stage.mat)



  return(list(delta=delta,covariance.delta=covariance.delta,second.stage.mat = second.stage.mat,second.stage.test,global.test,first.stage.mat,first.stage.test,loglikelihood = loglikelihood,
              AIC = AIC,beta=beta,covariance.beta=covariance.beta,
              z.standard=z.standard))

}
