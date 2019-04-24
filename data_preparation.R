
#packages + data upload


library(otrimle)

library(rpart)
library(monmlp)
library(kernlab)
#libraries
library(tseries)
library("tsensembler")

library("forecast")


library(kernlab)


library(factoextra)

library(NbClust)

set.seed(123)


#mars
library("mda")


### gbm
library(gbm)

###rf
library(randomForest)


library(RSNNS)


#gp gaussian processes
library(kernlab)


#pls
## add ncomp to predict
library(pls)


#####

#CONSTRUCT A SLIDING WINDOW ENSEMBLE


data=data('ice.river')


plot(prec)


fj=data.frame(flow.vat)

#data emnedding preparation


dataset=data_lag_prep(fj,1,5)



train.ts=data.frame(dataset[1:815,1])
test.ts=data.frame(dataset[816:1091,1])

names(train.ts)=names(test.ts)="target"

train.reg=dataset[1:815,]
test.reg=dataset[816:1091,]



