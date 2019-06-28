

path= '/Users/saadalla/Desktop/'

setwd(path)


data=read.csv('bike_sharing_registred_counts_ts.csv')


#data preparation and time series embedding preparation

dataset=data_lag_prep(data,1,5)


#split ts into train and test for regression and time series analysis models

train.ts=data.frame(dataset[1:2000,1])
test.ts=data.frame(dataset[2001:2688,1])

names(train.ts)=names(test.ts)="target"

train.reg=dataset[1:2000,]
test.reg=dataset[2001:2688,]


