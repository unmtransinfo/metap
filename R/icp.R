#############################################################
### ICP: Inductive Conformal Prediction for Classification
#        using Random Forest
#############################################################

#helper function
library(randomForest)



processFrameIntoMatrix = function(dataFrameInput) {

  #print(head(colnames(dataFrameInput),10))
  #stop()

  X <- data.matrix(dataFrameInput[,-c("Y", "id1", "subset")])  
  y <- ifelse(dataFrameInput[, Y] == "pos", 1, 0)
  dtrain <- xgb.DMatrix(X, label = y)
  return(dtrain)

}

processIdsIntoMatrix = function(data,negativeIds,positiveIds) {

  X <- data.matrix(data[id1 %in% c(positiveIds, negativeIds), -c("Y", "id1", "subset")])  
  y <- ifelse(data[id1 %in% c(positiveIds, negativeIds), Y] == "pos", 1, 0)
  dtrain <- xgb.DMatrix(X, label = y)
  return(dtrain)
}


CPObsFuzziness = function(matPValues, testLabels)
{
  if (is.null(matPValues) || is.null(testLabels)){
    stop("\n 'matPValues' and 'testLabels' are required as input\n")
  }

  nrTestCases = length(testLabels)

  sumPValues = 0
  for(indxTestSet in 1:nrTestCases)
  {
    exclude = testLabels[indxTestSet] #exclude the p-value of the true label
    sumPValues = sumPValues + sum(matPValues[indxTestSet, -exclude])
  }
  result = sumPValues/nrTestCases
  return(result)
}


#' Computes the deviation from exact validity as the Euclidean norm of
#' the difference of the observed error and the expected error
#' @param matPValues Matrix of p-values
#' @param testLabels True labels for the test-set
#' @return The deviation from exact validity
#' @export
CPValidity = function(matPValues = NULL, testLabels = NULL)
{
  if (is.null(matPValues) || is.null(testLabels)){
    stop("\n 'matPValues' and 'testLabels' are required as input\n")
  }

  signifSet = seq(.01, .99, by=.01) #significance level set

  nrTestCases = length(testLabels)
  errAtSignif = rep(0, length(signifSet))

  for(indx in 1: length(signifSet)){
    signifTest = (matPValues  > signifSet[indx])*1

    err = 0
    for(i in 1:nrTestCases)
    {
      err = err + ( (signifTest[i, testLabels[i]] == 0) * 1 )
    }
    err = err/nrTestCases
    errAtSignif[indx] = (err - signifSet[indx])^2
  }

  result = sqrt(sum(errAtSignif))
  return(result)
}


CPEfficiency = function(matPValues, testLabels, sigfLevel = 0.05)
{
  if (is.null(matPValues) || is.null(testLabels)){
    stop("\n 'matPValues' and 'testLabels' are required as input\n")
  }

  nrTestCases = length(testLabels) #size of the test set

  signifTest = (matPValues  > sigfLevel)*1 #compute the prediction region

  err = 0
  for(i in 1:nrTestCases)
  {
    err = err + ( (sum(signifTest[i, ]) > 1) * 1 )
  }
  

  print("DOUBLE VALUES")
  print(err)
  
  result = err/nrTestCases

  return(result)
}
CPErrorRate = function(matPValues, testLabels, sigfLevel = 0.05,RT)
{
  if (is.null(matPValues) || is.null(testLabels)){
    stop("\n 'matPValues' and 'testLabels' are required as input\n")
  }

  nrTestCases = length(testLabels)

  signifTest = (matPValues  > sigfLevel)*1  # if a value in the matrix is bigger than sig, set it to 1

  err = 0
  for(i in 1:nrTestCases)
  {

    booleanFlag = 2
    if(testLabels[i]=="pos") {
      booleanFlag= 1
    }


    err = err + (signifTest[i,booleanFlag] == 0)*1
  }

  result = err/nrTestCases
  return(result)

}


# fitModel = function(trainingSet=NULL, method = "rf",  nrTrees = 100)
# {
#   if(is.null(trainingSet))
#   {
#     stop("Error: training set is NULL")
#     return(NULL)
#   }

#   if(method != "rf")
#   {
#     stop("Error: only random forest is supported in the current release")
#     return(NULL)
#   }

#   #the first colum should be class labels, labeled as 1, 2, ...

#   names(trainingSet)[1] <- "Class"

#   #print(is.data.frame(data.frame(trainingSet)))
#   #print(is.data.frame(trainingSet))

#   trainingSet$Class <- as.factor(trainingSet$Class)
  
#   #trainingSet$X <- as.numeric(trainingSet$X)
#   #print(trainingSet)
#   #norm.votes = TRUE, keep.forest = TRUE, predict.all = TRUE

#   rfModel <- randomForest(Class ~ ex, data = trainingSet
#                           , type = "classification")

#   return(rfModel)
# }

computeConformityScores = function(modelFit = NULL, calibrationSet = NULL)
{

  #The first colum should be the class labels

  #calibLabels = as.numeric(as.factor(calibrationSet[, 1]))
  
  calibLabels = (as.numeric(unlist(calibrationSet[,1]))-1)

  #calibLabels = head(calibLabels,250)
  # calibration set is still dataframe here

  predProb = predict(modelFit, processFrameIntoMatrix(calibrationSet), type="prob") # run predictions on data matrix
  # this has got to come out correctly ....... 
  
  # we need to solve the score now, to get a computed / per class option 
  #predProb = data.table("0"=(1-predProb),"1"=predProb)

  positiveProb = predProb
  negativeProb = (1-predProb)
  nrLabels = ncol(predProb)-1 # number of class labels . here we will encode this differently 


  # one of the labels, w/ class members 
  #classMembers = which(calibLabels == i)

  positiveMembers = which(calibLabels == 1) # for each of these indexes, what the values in positivePROB

  negativeMembers = which(calibLabels == 0) # for each of these indexes, what are the values in negativePROB

  MCListConfScores = list(positiveProb[positiveMembers],negativeProb[negativeMembers]) #Moderian Class wise List of conformity scores
  
  #print(MCListConfScores[[1]])

  return(MCListConfScores)
}

computePValues = function(MCListConfScores = NULL, testConfScores = NULL)
{

  # given a list of the class values for each class MCLST
  # and an array of the predictions from model, modeled as dataframe where classes are on X axis

  #testConfScores = data.table("0"=(1-testConfScores),"1"=testConfScores)
  
  negativeModelScores = (1-testConfScores)
  positiveModelScores = testConfScores
  
  nrTestCases = length(testConfScores) # what is this thing now?? # this is right from model,


  pValues = matrix(0, nrTestCases,  2) # matrix created with the binary labels

  # for each test case, for each label, we compute the number of values which are <= to the corresponding score from the model
  # divide this, by the number of entries in the give class ... assign to that case, value= represents the probability of this value 

  for(k in 1:nrTestCases)
  {
      
    #class conf, is the score for that class ... 
    # alpha is the core our test conf got for this particular case 
    # pval is the length calculation 

    #alpha = testConfScores[k, ..l] # from model
    #pVal = length(which(classConfScores < alpha)) + runif(1) * length(which(classConfScores == alpha))

    aPos = positiveModelScores[k]
    PclassConfScores = MCListConfScores[[1]]
    pVal = length(which(PclassConfScores < aPos)) + runif(1) * length(which(PclassConfScores == aPos)) #runif(1) * 
    pValues[k,1] = pVal/(length(PclassConfScores)+1) # why add the plus 1 here?

    nPos = negativeModelScores[k]
    NclassConfScores = MCListConfScores[[2]]
    npVal = length(which(NclassConfScores < nPos)) + runif(1) * length(which(NclassConfScores == nPos)) #runif(1) * 
    pValues[k,2] = npVal/(length(NclassConfScores)+1) # why add the plus 1 here?

  }
  return(pValues)
}

#' Class-conditional Inductive conformal classifier for multi-class problems
#' @param trainingSet Training set
#' @param testSet Test set
#' @param ratioTrain The ratio for proper training set
#' @param method Method for modeling
#' @param nrTrees Number of trees for RF
#' @return The p-values
#' @export
# ICPClassification = function(trainingSet, testSet,  ratioTrain = 0.7, method = "rf", nrTrees = 100)
# {
#   if(is.null(trainingSet) || is.null(testSet) )
#   {
#     stop("\n 'trainingSet' and 'testSet' are required as input\n")
#   }

#   nTrainSize = nrow(trainingSet)
#   #create partition for proper-training set and calibration set.
#   result = sample(1:nTrainSize,  ratioTrain*nTrainSize)

#   result = c(1,2,3,4) 
  
#   #print(result) 
#   #print(typeof(result))

#   calibSet = trainingSet[-result, ]
  
#   properTrainSet = trainingSet[result, ]

#   modelFit = fitModel(properTrainSet, method = method,  nrTrees = nrTrees)
  
#   if(is.null(modelFit))
#     return(NULL)

#   MCListConfScores = computeConformityScores(modelFit, calibSet)

#   testConfScores = predict(modelFit, testSet[, -1], type = "prob")
  
#   print("===TEST_CONF===")
#   print(testConfScores)

#   pValues = computePValues(MCListConfScores, testConfScores)
#   print(pValues)
#   print("=====")
#   print("=====")
#   #print(pValues)

#   TestSmall <- data.frame("ex" = c(1),"g" = c(1)) # uses first val as label
#   #print(predict(modelFit, TestSmall, type="prob"))


#   predictions = predict(modelFit, testSet[, -1])

#   #print(predictions,testSet[,1])

#   return(pValues)
# }


# exData <- sample(1:10,140,replace=T) # create, 1000 values in range # if value is 1-3> true, if its 3 -> 8 false, if its greater than 8, true

# y <- c()
# index = 1
# for(n in exData) {
  
#   exData[index] = exData[index] * runif(1)
  
#   if(n < 1) {
#     y <- c(y,2)
#   } else {

#     randomNoise = runif(1)
#     if(randomNoise > 0.7) {
#       y <- c(y,2)
#       } else {
#       if(n > 9) {
#         y <- c(y,2)
#       } else {
#         y <- c(y,1)
#       }
#     }
#   }
#   index = index + 1
# }

# #print(y)
# #DATA <- data.frame("Y" = c(2,2,2,1,2,1,2,1,2,1,2 ,1,1,2,2 ),"ex" = c(1.3,2.11,2,1.9,2,1,2,1,2,1,2 ,1,1,2,2),"g" = c(1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1)) # uses first val as label
# #data 
# #print(colnames(DATA))

# DATA <- data.frame("Y" = y,"ex" = exData,"g" = sample(1:10,140,replace=T))

# train = head(DATA,70)
# test = tail(DATA,70) # bigger test

# pData <- ICPClassification(train,test,nrTrees = 0)

# #print(a[0])
# #modelItself <- a[1]
# #pData <- a[0]

# print("ERROR RATE, 70% prediction region")
# CPErrorRate(pData,test$Y,sigfLevel=0.30)

# print("EFFICIENCY RATE, 70% prediction region")
# CPEfficiency(pData,test$Y,sigfLevel=0.30)

# print("FUZZ,")
# CPObsFuzziness(pData,test$Y)

# print("DEVIATION FROM VALIDITY,")
# CPValidity(pData,test$Y)



#VALIDITY - errors happen less than 1 - E????? errors happen less than 95% ? so all the time..

#pValues, testSet, color="blue"

#sigLevel = 0.20
#for(row in 1:nrow(a)) {
#  for(cell in 1:ncol(a)) {  
#    cellVal = a[row,cell]
    
#    if(cellVal > sigLevel) {
#       a[row,cell] = "YES"
#    }
    #print(cellVal)
#  }
#}

