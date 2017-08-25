load("./data/simulation_ER_PR_HER2/simulated_data.Rdata")
y = data[[1]]
x = data[[2]]
colnames(y) <- c("Behavior","ER","PR","HER2","Grade")
HeterResultFunction <- function(y,
                                  baselineonly=NULL,
                                  main.effect=NULL,
                                  pairwise.interaction=NULL,
                                  saturated=NULL,
                                missingTumorIndicator = 888){
  tumor.number <- ncol(y)-1
  y.case.control <- y[,1]
  y.tumor <- y[,2:(tumor.numer+1)]
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
  z.design.main.effect <- GenerateZDesignMainEffect(tumor.character.cat,
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
  z.all <- ZDesigntoZall(baselineonly,
                         main.ffect,
                         pairwise.interaction,
                         saturated,
                         z.design.baselineonly,
                         z.design.main.effect,
                         z.design.pairwise.interaction,
                         z.design.saturated)
    delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
    x.all <- GenerateXAll(y,baselineonly,main.effect,pairwise.interaction,saturated)
    ###z standard matrix means the additive model z design matrix without baseline effect
    ###z standard matrix is used to match the missing tumor characteristics to the complete subtypes

    z.standard <- z.design.main.effect[,-1]
    y_em <- ProbFitting(delta0,y,x.all,z.standard,z.all,missingTumorIndicator)





}

source("./R/support_temp.R")

###Generate all the x

GenerateXAll <- function(y,baselineonly,main.effect,pairwise.interaction,saturated){
  n = nrow(y)
  x.all = matrix(1,n)
  if(is.null(baselineonly)==0){
    x.all = cbind(x.all,baselineonly)
  }else if(is.null(main.effect)==0){
    x.all = cbind(x.all,main.effect)
  }else if(is.null(pairwise.interaction)==0){
    x.all = cbind(x.all,pairwise.interaction)
  }else if(is.null(saturated)==0){
    x.all = cbind(x.all,saturated)
  }
  return(x.all)
}





####Calculate the conditional probability for E step
ProbFitting <- function(delta0,y,x.all,z.standard,z.all,missingTumorIndicator){
  result <- matrix(0,nrow=nrow(y),ncol = nrow(z.standard))
  beta <- matrix(z.all%*%delta0,ncol = nrow(z.standard))
  #beta <- matrix(z_all%*%delta0,nrow = ncol(z_design))
  #beta <- matrix(z_all%*%delta0,nrow = nrow(z_standard))
  for(i in 1:nrow(y)){
    if(y[i,1]==1) {
      idx <- which(y[i,2:ncol(y)]!=missingTumorIndicator)
      if(length(idx)==0) {
        result[i,] <- exp(x.all[i,]%*%beta)/(1+sum(exp(x.all[i,]%*%beta)))
      }else{
        jdx <- apply(z.standard,1,function(t){all(t[idx]==y[i,idx+1])})
        jdx <- which(jdx==T)
        temp <- exp(x.all[i,]%*%beta[,jdx])
        result[i,jdx] <- temp/sum(temp)
      }
    }
  }

  return(result)
}


