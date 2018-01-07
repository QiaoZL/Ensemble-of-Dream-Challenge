
setwd('D:/COURSE/bioinformatic/submissions')
options(warn =-1)

library(AUC)
library(ROCR)
library(pROC)

# get 5-fold results of each model

r1 <- read.csv("ViralChallenge_Time24_CLINICAL_round_1.csv",  header=T, row.names=1)
r2 <- read.csv("ViralChallenge_Time24_CLINICAL_round_2.csv",  header=T, row.names=1)
r3 <- read.csv("ViralChallenge_Time24_CLINICAL_round_3.csv",  header=T, row.names=1)
r4 <- read.csv("ViralChallenge_Time24_CLINICAL_round_4.csv",  header=T, row.names=1)
r5 <- read.csv("ViralChallenge_Time24_CLINICAL_round_5.csv",  header=T, row.names=1)
Lei_result.maintable <- list(r1,r2,r3,r4,r5)
#read.table("gene_expression_n438x978.txt", sep="\t", header=T, row.names=1)

#pre-poccessing
#parts 1, lei's submission

Lei_result.prefunc <- function(roundres,firstcolumn)
{
  result <-  data.frame(row.names = rownames(roundres),
                        possibility = vector(length = nrow(roundres)),
                        label = roundres[,1])
  
  for(i in 1:nrow(roundres))
  {
    numresult<-0
    sumresult<-0
    for(j in 0:4)
    {
      iter = firstcolumn+j
      if(!is.na(roundres[i,iter]))
      {
        sumresult = sumresult+roundres[i,iter]
        numresult = numresult+1
      }
    }
    result[i,'possibility'] <- sumresult/numresult
  }
  
  return (result)
  
}

mergefunc <- function(preres)
{
  finaltable = preres[[1]]
  for(i in 2:length(preres))
  {
    preres[[i]]$label = NULL
    finaltable <- merge(finaltable, preres[[i]],
                                 by = 'row.names')
    finaltable <- transform(finaltable,row.names=Row.names, Row.names=NULL)
  }
  
  finaltable$finalP <- rowMeans(finaltable[,!(names(finaltable)=='label')])
  finaltable <- data.frame(finaltable[,colnames(finaltable) %in% c('label','finalP')])
    
  return(finaltable)
}

tablemergefunc <- function (modellist)
{
  print(names(modellist[[1]]))
  names(modellist[[1]])[2] = "model 1"
  finaltable = modellist[[1]]
  
  for(i in 2:length(modellist))
  {
    modellist[[i]]$label = NULL
    print(names(modellist[[i]]))
    names(modellist[[i]])[1] <- paste('model',i)
    finaltable <- merge(finaltable, modellist[[i]],
                        by = 'row.names')
    finaltable <- transform(finaltable,row.names=Row.names, Row.names=NULL)
  }
  
  #finaltable$finalP <- rowMeans(finaltable[,!(names(finaltable)=='label')])
  #finaltable <- data.frame(finaltable[,colnames(finaltable) %in% c('label','finalP')])
  
  return(finaltable)
}

firstcolumn<-c(2,7,12,17)
Lei_result.pocessed.lr<-list()
Lei_result.pocessed.svm<-list()
Lei_result.pocessed.knn<-list()
Lei_result.pocessed.mnb<-list()
for (i in 1:5)
{
  Lei_result.pocessed.lr<-c(Lei_result.pocessed.lr,list(Lei_result.prefunc(Lei_result.maintable[[i]],2)))
  Lei_result.pocessed.svm<-c(Lei_result.pocessed.svm,list(Lei_result.prefunc(Lei_result.maintable[[i]],7)))
  Lei_result.pocessed.knn<-c(Lei_result.pocessed.knn,list(Lei_result.prefunc(Lei_result.maintable[[i]],12)))
  Lei_result.pocessed.mnb<-c(Lei_result.pocessed.mnb,list(Lei_result.prefunc(Lei_result.maintable[[i]],17)))
}

Lei_result.final.lr <- mergefunc(Lei_result.pocessed.lr)
Lei_result.final.svm <- mergefunc(Lei_result.pocessed.svm)
Lei_result.final.knn <- mergefunc(Lei_result.pocessed.knn)
Lei_result.final.mnb <- mergefunc(Lei_result.pocessed.mnb)

#parts2 Li's submission

Li_result.pocessed.lasso <- list()
Li_result.pocessed.rm <- list()
Li_result.pocessed.lassosvm <-list()
Li_result.pocessed.genesetsvm <- list()

unifycolnames <- c('possibility', 'label')
for (i in 1:5)
{
  filename1 <- paste('result_LOOCV_LASSO_glm_phase1',i,'.csv')
  filename2 <- paste('result_LOOCV_GeneSet_glm_phase1',i,'.csv')
  filename3 <- paste('result_LOOCV_LASSO_SVM_phase1',i,'.csv')
  filename4 <- paste('result_LOOCV_GeneSet_SVM_phase1',i,'.csv')
  
  Li_result.pocessed.lasso <-c( Li_result.pocessed.lasso ,list(read.csv(filename1,  header=T, row.names=1)))
  colnames(Li_result.pocessed.lasso[[i]]) <- unifycolnames
  
  Li_result.pocessed.rm <-c( Li_result.pocessed.rm ,list(read.csv(filename2,  header=T, row.names=1)))
  colnames(Li_result.pocessed.rm[[i]]) <- unifycolnames
  
  Li_result.pocessed.lassosvm <-c( Li_result.pocessed.lassosvm ,list(read.csv(filename3,  header=T, row.names=1)))
  colnames(Li_result.pocessed.lassosvm[[i]]) <- unifycolnames
  
  Li_result.pocessed.genesetsvm <-c( Li_result.pocessed.genesetsvm ,list(read.csv(filename4,  header=T, row.names=1)))
  colnames(Li_result.pocessed.genesetsvm[[i]]) <- unifycolnames
  
}

Li_ressult.final.lasso <- mergefunc(Li_result.pocessed.lasso)
Li_ressult.final.rm <- mergefunc(Li_result.pocessed.rm)
Li_ressult.final.lassosvm <- mergefunc( Li_result.pocessed.lassosvm)
Li_ressult.final.genesetsvm <- mergefunc(Li_result.pocessed.genesetsvm)

table.mergelist <- list(Li_ressult.final.lasso,Li_ressult.final.rm,Li_ressult.final.lassosvm,Li_ressult.final.genesetsvm,
                        Lei_result.final.lr,Lei_result.final.svm,Lei_result.final.knn,Lei_result.final.mnb)
table.final <- tablemergefunc (table.mergelist)

#build ensemble model
ensemble.means <- table.final
ensemble.means$means <- rowMeans(ensemble.means[,2:ncol(ensemble.means)])

ensemble.vote <- table.final
ensemble.vote$median <- rowMedians(as.matrix(ensemble.vote[,2:ncol(ensemble.vote)]))

ensemble.glm.traindata <- table.final
ensemble.glm.model <- glm(label~.,data = ensemble.glm.traindata)
ensemble.glm.traindata$predict_result = predict(ensemble.glm.model ,
                                                ensemble.glm.traindata[,2:ncol(ensemble.glm.traindata)], 
                                                type='response')





