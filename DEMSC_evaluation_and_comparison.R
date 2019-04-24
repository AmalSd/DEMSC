



#example of evaluation

pred.a=calculate_arima_prediction(train.ts, test.ts,nrow(test.reg),1)
rmse(pred.a,test.reg$target)
plot(test.reg$target,type='l')
lines(pred.a,col='red')



data.train=train.reg
data.test=test.reg
n=nrow(data.train)
val.length=50
data.val=data.train[((n-val.length+1):n), ]
data.train.n=data.train[(1:(n-val.length)), ]
formula=target~.
models=train.models(data.train.n, formula,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))

predictions.table=predict.models(models,data.test,formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))


sapply(1:ncol(predictions.table), function(x) rmse(predictions.table[,1],predictions.table[,x]))


tp1=15
lim=0.1
val.length=20

updated_selection=topk.model.sel(models,data.train, data.test,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),lim,tp1)

updated_selection1=updated_selection$models.sel

alarm=updated_selection$alarm
H=8


pred.top=ens.top.pred(data.test,predictions.table,updated_selection1,H)



pred.ens.all.sw=sapply(1:nrow(data.test),function(t)ens.step.all(data.test,predictions.table,t,H))

pred.ens.all=sapply(1:nrow(data.test),function(t)ens.avg.all(data.test,predictions.table,t))



rmse(data.test$target,pred.top)
rmse(data.test$target,pred.ens.all)

rmse(data.test$target,pred.ens.all.sw)


val.length1=10
cluster.res=update.cluster(models,data.train, data.test,val.length1, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection,tp1)

pos.cl=lapply(1:length(cluster.res) ,function(x) if(length(cluster.res[[x]]$list.position)>=4){
  cluster.res[[x]]$list.position[1:4]
}else{c(cluster.res[[x]]$list.position[1],rep(cluster.res[[x]]$list.position[2],3))})


pred.cl=ens.top.pred(data.test,predictions.table,pos.cl,H)

st=Sys.time()
val.length1=10
tp.cl.all=compute.cluster.tp.all(models,data.train,data.test,val.length1,H, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection)
ed=Sys.time()

pred.tp.cl.all=tp.cl.all[[1]]

pos.models=tp.cl.all[[2]]



val.length1=10

pred.cl.or=compute.cluster.or.all(models,data.train, data.test,val.length1,H, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))





val.length1=10
pred.tw.tp=compute.cluster.tw.tp.all(models,data.train, data.test,t,val.length1, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection)



st=Sys.time()
val.length1=50
pred.st.all=compute.cl.tp.stac.all(models,data.train, data.test,t,val.length1, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),pos.models,updated_selection)
ed=Sys.time()


sapply(1:12, function(x) rmse(data.test$target,pred.st.all[[x]]))

pred.cl.st=pred.st.all[[10]]





#stacking framework using Rf

t=1

val.length=400
output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))


sel.models=(output.m)

x= sel.models[,2:ncol(sel.models)]
y= sel.models[,1]
train.st1=sel.models
formula1=ts~.

test.st=predictions.table
names(test.st)=names(train.st)
x.test= test.st[,2:ncol(sel.models)]
y.test= test.st[,1]

fitControl <- trainControl(
  method = 'cv',                   # k-fold cross validation
  number = 5,                      # number of folds
  savePredictions = 'final'      # saves predictions for optimal tuning parameter
  
) 

mlboost=train(formula1,train.st1,method='rf',  trControl = fitControl)
pred.stacking=predict(mlboost,test.st,ncomp=1)
rmse(pred.stacking,data.test$target)


# Opera framework 
library(opera)
t=1

val.length=400


output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))

X=as.matrix(output.m[,2:ncol(output.m)])

Y=as.vector(output.m[,1])

X.new=as.matrix(predictions.table[,2:ncol(predictions.table)])

Y.new=as.vector(predictions.table[,1])

MLpol <- mixture(Y = Y, experts = X, model = "MLpol", loss.type = "square")


pred.mlpol=predict(MLpol,newexperts=X.new,newY=Y.new, type='response')

rmse(pred.mlpol,data.test$target)
fs <- mixture(Y = Y, experts = X, model = "FS", loss.type = "square")

pred.fs=predict(fs,newexperts=X.new,newY=Y.new, type='response')

rmse(pred.fs,data.test$target)

rmse(pred.tp.cl.all,data.test$target)


ewa <- mixture(Y = Y, experts = X, model = "EWA", loss.type = "square")

pred.ewa=predict(ewa,newexperts=X.new,newY=Y.new, type='response')

rmse(pred.ewa,data.test$target)

ogd <- mixture(Y = Y, experts = X, model = "OGD", loss.type = "square")

pred.ogd=pred.ewa=predict(ogd,newexperts=X.new,newY=Y.new, type='response')


rmse(pred.ogd,data.test$target)




#Ade framework cmp

st1=Sys.time()

library(tsensembler)

specs=model_specs(learner = c("bm_ppr", "bm_svr", "bm_randomforest",
                              "bm_gaussianprocess",
                              "bm_gbm", "bm_pls_pcr", "bm_ffnn", "bm_mars"),
                  learner_pars = list( bm_ppr = list(nterms = c(2,4),
                                                     sm.method = "supsmu" ),
                                       bm_svr = list(kernel =  c("vanilladot", "polydot","rbfdot"),
                                                     C = c(1),epsilon = .01),
                                       bm_randomforest = list(
                                         num.trees = c(500,250,100) ),
                                       bm_gbm = list(
                                         interaction.depth = 1,
                                         shrinkage = c(.01),
                                         n.trees = c(5,10)
                                       ),
                                       bm_mars = list(
                                         nk = c(7,15),
                                         degree = c(3,1),
                                         thresh = .001
                                       ),
                                       bm_ffnn = list(
                                         size = c(3,5,7,10,15,25),
                                         decay = .01
                                       ),
                                       bm_pls_pcr = list(
                                         method = c("kernelpls", "simpls")
                                       ),
                                       bm_gaussianprocess = list(
                                         kernel =  c("vanilladot", "polydot","rbfdot"),
                                         tol = .01)))


ade <- ADE(target ~., data.train, specs)

pred.ade=predict(ade,data.test)

rmse(pred.ade@y_hat,data.test$target)

ed1=Sys.time()


rmse(pred.stacking,data.test$target)
dets <- DETS(target ~., data.train, specs, lambda = 30, omega = .2)

pred.dets=predict(dets,data.test)

rmse(pred.dets@y_hat,data.test$target)


