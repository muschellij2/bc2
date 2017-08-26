load("simulated_data.Rdata")
y = simulated_data[[1]]
x = simulated_data[[2]]
y.pheno.complete <- simulated_data[[3]]
x.pheno.complete <- simulated_data[[4]]
colnames(y.pheno.complete)<- colnames(y) <- c("Behavior","ER","PR","HER2","Grade")
colnames(x.pheno.complete) <- colnames(x) <- c("PC1","PC2")
try = matrix(rnorm(nrow(y),0),nrow(y),1)
missingTumorIndicator = 888







library(devtools)
install_github("andrewhaoyu/bc2")
library(bc2)


result4 = EMmvpoly(y,baselineonly = NULL,main.effect = x,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)










































source("./mvpoly_no_missing.R")
result.nomissing <- MvpolyNoMissing(y.pheno.complete,
                                    baselineonly=NULL,
                                    main.effect=x.pheno.complete,
                                    pairwise.interaction=NULL,
                                    saturated=NULL)







colnames(try) <- c("gene")
write.csv(try,"gene.csv",row.names=F,quote=F)

y <- read.csv("y.csv",header=T)
x <- read.csv('x.csv',header=T)
gene <- read.csv("gene.csv",header=T)
source('~/GoogleDrive/breast_cancer/V10/known_SNPs_anlysis/try/try.R')
result = HeterResultFunction(y,baselineonly = try,main.effect = x,pairwise.interaction = NULL,saturated = NULL)
source('./try2.R')
result2 = HeterResultFunction(y,baselineonly = try,main.effect = x,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
source('./try3.R')
result3 = HeterResultFunction(y,baselineonly = try,main.effect = x,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
source('./try4.R')
result4 = HeterResultFunction(y,baselineonly = gene,main.effect = x,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

source('./try5.R')
result5 = ScoreTestSupport(y,baselineonly = try,main.effect = x,pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

source('./try6.R')
result6 = ScoreTest(y,try,second.stage.structure = "main.effect",score.test.support=result5,missingTumorIndicator=888)
















library(devtools)
install_github("andrewhaoyu/bc2")
library(bc2)
library(readr)





y_standard <- as.matrix(read_csv("./onco_y_standard.csv"))
x <- as.matrix(read_csv("./gene_and_pc.csv"))
y_onco <- as.matrix(read_csv("./onco_y.csv"))
pc <- x[,2:11]
gene <- x[,1]


##############Two stage model MLE
Heter_result = EMmvpoly(y_standard,baselineonly = NULL,main.effect = as.matrix(cbind(x)),pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)


Heter_result_G = EMmvpoly(y_onco,baselineonly = NULL,main.effect = as.matrix(cbind(x)),pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)


#############Score Support Function
Score_support = ScoreTestSupport(y_onco,baselineonly = as.matrix(gene),main.effect = as.matrix(pc),pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)


############Score Test

Score.test.result = ScoreTest(y_onco,gene,second.stage.structure = "main.effect",score.test.support=Score_support,missingTumorIndicator=888)




