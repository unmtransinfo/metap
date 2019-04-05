library(R6)

library(data.table)
library(xgboost)
library(Matrix)
library(pROC)
library(ROCR)
library(PRROC)

# Person <- R6Class(classname="Context",public = list(
# 	name = NULL,
# 	age = NULL,
# 	initialize = function(name=NA,age=NA) {
# 		self$name = name
# 		self$age = age
# 	},
# 	greet = function() {
# 		cat(paste0("hello my name is ",self$name,"\n"))
# 	}
# ))

Output <- R6Class(classname="Output",public = list(
	DataObject = NULL,
	initialize = function(Data) {
		self$DataObject = Data
	},
	AUCROC = function() {
		data = self$DataObject
		frame = data$d
		# this will be recorded... with the hyper parameters?
		#print()
		print(auc(as.numeric(frame[,frame$label]),frame[,frame$pred]))
	},
	AUCPR = function() {

		data = self$DataObject
		frame = data$d
		plot(pr.curve(frame[,frame$pred],as.numeric(frame[,frame$label]),curve=TRUE))
	}
))
Context <- R6Class(classname="Context",public = list(
	experiment = "test",
	initialize = function() {
		
	},
	run = function() {
		cat(paste0("hello my name is ",self$experiment,"\n"))
	}
))


Data <- R6Class(classname="Data",public = list(
	d = NULL,
	labelString = NULL,
	hasLabel = TRUE,
	possibleLabels = NULL,
	initialize = function(dataFrame=NULL,label=NULL) {
		if(!is.null(dataFrame) & !is.null(label)) {
			self$loadAsDataFrame(dataFrame)
			self$makeLabel(label)
		}
	},
	loadAsDataFrame = function(df) {
		self$d = df
	},
	loadFromFile = function(type,path) {
		# switch on type
		RDS = readRDS(path)
		self$loadAsDataFrame(RDS)
	},
	only = function(label) {
		# you can still split it here...
		newD = self$d[Y == label]
		newSet <- Data$new(newD,self$labelString)
		return(newSet)
	},

	# dropRowsWhere ::

	removeAllBut = function(...) {
		keep <- c(...)
		self$d = self$d[,..keep]
	},
	removeFeatures = function(...) {
		drop <- c(...)
		self$d = self$d[,-..drop]
	},
	countFeatures = function() {
		return (ncol(self$d))
	},
	makeLabel = function(labelString) {
		self$labelString = labelString
		self$possibleLabels = unique(self$getLabelData())
	},
	getLabelData = function() {
		return (self$d[,get(self$labelString)])
	},
	unLabeledSet = function() {
		newSet <- Data$new(self$d,self$labelString)
		newSet$hasLabel = FALSE
		newSet$removeFeatures(self$labelString)
		return(newSet)
	},
	count = function() {
		return (nrow(self$d))
	},
	getDR = function() {
		return(self$d)
	},
	shuffle = function() {

	},
	sampleSplit = function(splitRatio,index=NULL) { # samples evenly from both classes... lets return a dictionary here ... which will contain our split 

		posIds = NULL

		copy = self$getDR()

		splitFrame = NULL
		remainingFrame = NULL

		originOrder = colnames(self$d)
		
		for(i in self$possibleLabels) {

			classSlice = copy[copy$Y == i]

			inTrain <- classSlice[,sample(.N, floor(.N*splitRatio))] # sample from this class slice
			
			Train <- classSlice[inTrain]
			Test <- classSlice[-inTrain]

			if(is.null(splitFrame)) {
				splitFrame = Train
			} else {
				splitFrame = rbind(splitFrame,Train)
			}

			if(is.null(remainingFrame)) {
				remainingFrame = Test
			} else {
				remainingFrame = rbind(remainingFrame,Test) #rbind(remainingFrame,Test)
			}
			
		}

		result = c()

		splitSet <- Data$new()
		remainSet <- Data$new()
		
		splitSet$loadAsDataFrame(splitFrame[,..originOrder])
		remainSet$loadAsDataFrame(remainingFrame[,..originOrder])
		
		splitSet$makeLabel(self$labelString)
		remainSet$makeLabel(self$labelString)

		result[["split"]] = splitSet
		result[["remaining"]] = remainSet


		return(result)

	}
))

# example APIs for all kinds of stuff ... GOOD version of abstract class
XGBModel <- R6Class(classname="XGBModel",public = list(
	model = NULL,
	initialize = function(trainData) {
		
		tr = self$preptrain(trainData)

		#print(head(colnames(tr),3))
		result = xgboost(
		    data=tr,
		    max_depth=10,
		    eta = 0.15,
		    gamma = 0.01,
		    min_child_weight = 0, 
		    subsample = 1,
		    colsample_bytree = 0.8,
		    objective = "binary:logistic",
		    nrounds = 2,
		    verbose = 2,
		    eval_metric='auc'
		)

		self$model = result

	},
	preptrain = function(train) { # functions for each model, models may have specifc reqs 

		noLabels = train$unLabeledSet()
		X <- data.matrix(noLabels$getDR())
		y <- ifelse(train$getLabelData() == "pos",1,0)
		return(xgb.DMatrix(X, label = y))

	},
	prepvalid = function() { # interface

	},
	predict = function(test) { # interface
		noLabels = test$unLabeledSet()
		d = xgb.DMatrix(data.matrix(noLabels$getDR()))
		pred <- predict(self$model,d)

		cl = Data$new()

		if(test$hasLabel) {
			labels = test$getLabelData()
			cl$loadAsDataFrame(data.table("pred"=pred,"label"=labels))
		} else {
			cl$loadAsDataFrame(data.table("pred"=pred))
			cl$hasLabel = FALSE
		}	

		return(cl)
	},
	crossValidate = function(startData) { # interface

	}
))

# OUTPUT IS ANOTHER KIND OF DATA WRAPPER, WITH FUNCTIONS LIKE PR CURVES AND SO FORTH....
# RESULTS MIGHT BE MORE LIKE FUNCTIONS ... explicitly
# y <- ifelse(dt[id1 %in% c(train.pos.idx, train.neg.idx), Y] == "pos", 1, 0)


#Experiment <- Context$new()
#Experiment$run()
#Context$run()


# we'd build up the logging system





# what do you want to implement first???

# like... divide data sets into 30% each 
# = 30% train / 30% split / 30% validation 








# ### SQL:  where | order by   select | update  group by

# #DF[with(DF, x > 1), ]

# stop()

# flights = fread("testdata.csv")

# print(is.data.table(flights))
# print(dim(flights))

# d = flights
# adep_delay = c("dep_delay","arr_delay")
# print(adep_delay)
# ans = flights[,-..adep_delay]
# head(ans)

# ## apparently this can be done??

# inTrain <- d[,sample(.N, floor(.N*.90))]
# Train <- d[inTrain]
# Test <- d[-inTrain]

# print(dim(Train))
# print(dim(Test))





#given a list of ids 


#print(rownames(flights))

#ans <- flights[, .(arr_delay, dep_delay)]

# we could recompute these splits X number of times to get the monte carlo behavior 
#trainData = split[["split"]]
#print(nrow(trainData))

#print(split[["split"]]$count())
#print(split[["remaining"]]$count())

#nolabels = trainData$getLabelData()
#print(nolabels)

#classifier <- XGBModel$new(trainData)

# print(nolabels$countFeatures())
# print(trainData$countFeatures())
# print(nolabels$hasLabel)
# print(trainData$hasLabel)

# train on the sample, eval on the real .... save the results 




# df <- data.frame(row.names = c("apple", "banana", "orange", "lemon", "lime"), 
#              value = c(1:5))
# remove_these <- c("apple", "orange")
# #Now we find the indicies of the rows that need to be removed

# rows_to_remove <- which(row.names(df) %in% remove_these)
# # And use the same technique you were trying to use before to remove rows.
# print(df)
# df <- df[-rows_to_remove,]

# print(df)