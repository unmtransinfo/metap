library(R6)
source("Object.R")

library(data.table)
library(xgboost)
library(Matrix)
library(pROC)
library(ROCR)


ALSet <- Data$new()
ALSet$loadFromFile("rds","../ALsmall.rds")
ALSet$makeLabel("Y")
ALSet$removeFeatures("subset","id1") #"id1"

# we could then walk over 1000 values and do this split, train a model and save... all in a couple lines of cod e
split = ALSet$sampleSplit(0.8) # this will give a split, and save the rest of the results so those can be accessed too 
train = split$split
test = split$remaining

classifier <- XGBModel$new(train) # will build and train the model right away
results = classifier$predict(test)
output = Output$new(results)
output$AUCPR()
output$AUCROC()

neg = train$only("neg") # we can select only a class here...

negSplit = neg$sampleSplit(0.1)
print(negSplit$split$count())






# now we'd add conformal prediction to this system, get the confidence intervals, 
# and then add a way to stack in hyperparameters ... 
