EMmvpolySelfDesign <- function(y,
                    x,
                    z.design,
                     missingTumorIndicator = 888,
                    z.all=NULL){
  if(is.null(z.all)){

  }
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
  z.design.additive <- GenerateZDesignMainEffect(tumor.character.cat,
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
  #   pxx = EM.result[[3]]
  #   y_em = EM.result[[4]]
  #  score_support_result <- score_support(pxx,x.all,baselineonly,z.all,z.standard,y_em)
  #  #return(score_support_result)
  # score_test_mis_result <- score_test_mis(y_em,baselineonly,score_support_result)

  return(EM.result)
  #return(list(score_c=score_test_mis$score_c,infor_c = score_test_mis$infor_c))
  #return(EM.result)

}
