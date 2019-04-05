#!/usr/bin/env Rscript
##
#library(RPostgreSQL)
library(data.table)
library(xgboost)
library(Matrix)
library(pROC)
library(ROCR)
library(PRROC)

source("icp.R")



# load the AL data, 

# train 100 models with folds, 

# output the average AUC, across them all 

# that's it - no grid search 





fitModel = function(trainPosIDS,trainNegIDS)
{

  
  #X <- data.matrix(dt[id1 %in% c(trainPosIDS, trainNegIDS), -c("Y", "id1", "subset")])  
  y <- ifelse(dt[id1 %in% c(trainPosIDS, trainNegIDS), Y] == "pos", 1, 0)

  #dtrain <- xgb.DMatrix(X, label = y)
  
  dtrain <- processIdsIntoMatrix(dt,trainNegIDS,trainPosIDS)
  
  sumpos <- sum(y == 1)
  sumneg <- sum(y == 0)

  result = xgboost(
    data=dtrain,
    max_depth=10,
    eta = 0.15, 
    gamma = 0.01,
    min_child_weight = 0, 
    subsample = 1, 
    colsample_bytree = 0.8,  
    objective = "binary:logistic",
    nrounds=50,
    verbose = 2,
    eval_metric='auc'
  )

  #print(result)

  return(result)
}

ICPClassification = function(dataSET,trainPosSet,trainNegSet, testSet,  ratioTrain = 0.7, method = "rf", nrTrees = 100)
{
 

  #nTrainSize = nrow(trainingSet)
  #create partition for proper-training set and calibration set.
  #result = sample(1:nTrainSize,  ratioTrain*nTrainSize)

  #calibSet = tail(dataSET)
  #trainingSet[-result, ] # tail from the second , passed in 
  # tail from the other 
  calibPos = tail(trainPosSet,20)
  calibNeg = tail(trainNegSet,1000)

  # use this calib set for the labels 
  #properTrainSet = trainingSet[result, ]

  trainPosSet <- trainPosSet[!trainPosSet %in% calibPos]
  trainNegSet <- trainNegSet[!trainNegSet %in% calibNeg] # this is proper train set 

  modelFit = fitModel(trainPosSet,trainNegSet)

  calibSet = dataSET[(id1 %in% c(calibPos, calibNeg)),-c("id1","subset")]

  # we need to turn the calibSet into... XGB data 
  #print(colnames(calibSet)) # works, Y is label, is first 
  # we build the calib set here, reduce the other colummns, and plug it in 
  
  if(is.null(modelFit))
    return(NULL)

  MCListConfScores = computeConformityScores(modelFit, calibSet)

  TM = processFrameIntoMatrix(testSet)
  # now we need to prepare the test set...
  testConfScores = predict(modelFit, TM, type = "prob")

  print("===TEST_CONF===")

  # MCLIST is cominng out wrong, needs to be sorted by class, 

  pValues = computePValues(MCListConfScores, testConfScores)
  
  #TestSmall <- data.frame("ex" = c(1),"g" = c(1)) # uses first val as label
  #print(predict(modelFit, TestSmall, type="prob"))
  #predictions = predict(modelFit, testSet[, -1])
  #print(predictions,testSet[,1])

  return(list(modelFit,pValues))

}



fn <- 'ALsmall.rds'

protein <- read.csv("proteinList.csv")
setDT(protein)

#subset only the training set
dt <- readRDS(fn) # kill everything but the good features here..

# here's a sample - not needed for XGB

#set.seed(100)

train.pos.idx <- sample(dt[Y == "pos", id1], size = 1.0*dt[Y == "pos", .N])
train.neg.idx <- sample(dt[Y == "neg", id1], size = 1.0*dt[Y == "neg", .N])

print("TRAINING ON POS")
print(length(train.pos.idx))

#testData = dt[!(id1 %in% c(train.pos.idx, train.neg.idx)),]
testData = dt


# we'd pass this into ICP 

# calib set, represents back third of the inputs...
# take the length of the sample, 

print("NROWS TEST")
print(nrow(testData))
CPOUT <- ICPClassification(dt,train.pos.idx,train.neg.idx,testData)

# once we've got the train / test, we need to pass it into XGBoost 


theModel = CPOUT[[1]]
pData = CPOUT[[2]]

labels = testData$Y
raw = data.matrix(testData[,-c("Y","id1","subset")]) # test data no labels 
pred <- predict(theModel,raw)

resultTable = data.frame("labels"=labels,"pred"=pred,"CPSCORE"=pData[,1],"CNSCORE"=pData[,2])

#print("PREDICTED POSITIVES ABOVE SIG LEVEL")
sig = 0.05
#print(nrow(resultTable[resultTable$CPSCORE > sig,]))
#print("PREDICTED NEGS ABOVE SIG LEVEL")
#print(nrow(resultTable[resultTable$CNSCORE > sig,]))

resultTable$CLASSVAL[resultTable$CPSCORE > sig & resultTable$CNSCORE > sig] <- "both"
resultTable$CLASSVAL[resultTable$CPSCORE > sig & resultTable$CNSCORE < sig] <- "pos"
resultTable$CLASSVAL[resultTable$CNSCORE > sig & resultTable$CPSCORE < sig] <- "neg"
resultTable$CLASSVAL[resultTable$CNSCORE < sig & resultTable$CPSCORE < sig] <- "empty"

#print("TOTAL TARGETS")
#print(nrow(resultTable[resultTable$labels == "pos",]))

print("NUMBER OF EMPTY ROWS::")
print(nrow(resultTable[resultTable$CLASSVAL == "empty",]))

print("MISTAKES - negatives as positive")
print(nrow(resultTable[resultTable$CLASSVAL == "neg" & resultTable$labels == "pos",]))

print("POSITIVES - positives as negatives")
print(nrow(resultTable[resultTable$CLASSVAL == "pos" & resultTable$labels == "neg",]))


#print("TARGEST IN THE EMPTY::")
#print(nrow(resultTable[resultTable$CLASSVAL == "empty" & resultTable$labels == "pos",]))




#print("NUMBER OF BOTH ROWS::")
#print(nrow(resultTable[resultTable$CLASSVAL == "both",]))

#print("TARGETS IN BOTH::")
#print(nrow(resultTable[resultTable$CLASSVAL == "both" & resultTable$labels == "pos",]))

#print("NUMBER OF STRONG POS::")
#print(nrow(resultTable[resultTable$CLASSVAL == "pos",]))

#print("TARGETS IN STRONG POS::")
#print(nrow(resultTable[resultTable$CLASSVAL == "pos" & resultTable$labels == "pos",]))

#print(" ================= ================= ================= ")
#print("NUMBER OF STRONG NEG::")
#print(nrow(resultTable[resultTable$CLASSVAL == "neg",]))


print("AUC")
print(auc(labels,pred))
pr <- pr.curve( labels, pred )
print(pr)

PRINT_SHELF <- 0

if(PRINT_SHELF == 1) {

  sortedTable = resultTable[order(-pred),]
  DT = sortedTable

  topshelf = head(sortedTable,20)
  print(topshelf[topshelf$labels=="pos",])
  print("XGB DRUGT IN SORTED SHELF BASIC SORT ON ORDER")
  print(nrow(topshelf[topshelf$labels=="pos",]))

  print("===== STARTING CP ========")
  print("===== STARTING CP ========")

  topshelf = resultTable

  topshelf$CLASSVAL = "empty"

  sig = 0.50

  topshelf$CLASSVAL[topshelf$CPSCORE > sig & topshelf$CNSCORE > sig] <- "both"
  topshelf$CLASSVAL[topshelf$CPSCORE > sig & topshelf$CNSCORE < sig] <- "pos"
  topshelf$CLASSVAL[topshelf$CNSCORE > sig & topshelf$CPSCORE < sig] <- "neg"

  print("THIS MANY POSITIVE LABELS")
  print(length(topshelf$CLASSVAL[topshelf$CLASSVAL=="pos"]))
  print("THIS MANY TARGETS IN POSITIVE SECTION")
  print(rownames(topshelf[topshelf$CLASSVAL=="pos" & topshelf$labels=="pos",]))

  print("DR RATIO")
  print(length(topshelf[topshelf$CLASSVAL=="pos" & topshelf$labels=="pos",])/length(topshelf$CLASSVAL[topshelf$CLASSVAL=="pos"]))

  print("THIS MANY NEGATIVE GUESSES")
  print(length(topshelf$CLASSVAL[topshelf$CLASSVAL=="neg"]))

  print("THIS TARGETS IN NEGATIVE SECTION")
  print(rownames(topshelf[topshelf$CLASSVAL=="neg" & topshelf$labels=="pos",]))

  print("THIS IS BOTH")
  print(length(topshelf$CLASSVAL[topshelf$CLASSVAL=="both"]))

  print("THIS TARGETS IN BOTH SECTION")
  print(rownames(topshelf[topshelf$CLASSVAL=="both" & topshelf$labels=="pos",]))

  print("THIS IS EMPTY")
  print(length(topshelf$CLASSVAL[topshelf$CLASSVAL=="empty"]))

  print("THIS TARGETS IN EMPTY SECTION")
  print(rownames(topshelf[topshelf$CLASSVAL=="empty" & topshelf$labels=="pos",]))

  sortedTable = resultTable[order(-resultTable$CPSCORE),]
  topshelf = head(sortedTable,20)
  print(topshelf)
  print("CP--DRUGT IN SORTED SHELF")
  print(nrow(topshelf[topshelf$labels=="pos",]))
  print("CP--TOTAL DRUGT")
  print(nrow(resultTable[resultTable$labels=="pos",]))

}
# Can we extract the both values .... 

print("===========")

print("===== EVAL METRICS LATER =====")


print("ERROR RATE, 95% prediction region")
CPErrorRate(pData,labels,sigfLevel=0.05)

print("EFFICIENCY RATE, 95% prediction region")
CPEfficiency(pData,labels,sigfLevel=0.05)

print("DEVIATION FROM VALIDITY,")
CPValidity(pData,labels)


print(resultTable[resultTable$labels=="pos",])
# warnings()
# print("FUZZ,")
# CPObsFuzziness(pData,labels)
