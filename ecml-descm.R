#Functions
data_lag_prep=function(data,col.indx,k1){  
  # Pre-processig and preparibg data
  reg_data=(t(sapply((k1+1):(dim(data)[1]),function(x) return(data[c(x,(x-k1):(x-1)),col.indx]))))
  colnames(reg_data)= c('target',paste("lag",k1:1,sep=''))
  reg_data=as.data.frame(reg_data)
  return(reg_data)
}

set.seed(2167)


calculaModeloArima<-function(timeseries,th=14*48)
{
  
  #if th<15 days, produce warning
  if (length(timeseries)<th)
  {
    print(length(timeseries))
    print(th)
    print("WARNING: The supplied series has a size smaller than expected !!!")
  }
  
  #Select the model
  fit <- tryCatch(auto.arima(timeseries,allowdrift=FALSE,seasonal = FALSE), error=function(e) e)
  if (is.list(fit))
  {
    arma<-fit$arma
    myOrder<-c(arma[1],arma[6],arma[2])
    if (is.vector(myOrder) && length(myOrder)==3)
    {
      if ((myOrder[1]==4 && myOrder[2]==0 && myOrder[3]==4) || (myOrder[1]==0 && myOrder[2]==0 && myOrder[3]==0) || (myOrder[1]==4 && myOrder[2]==0 && myOrder[3]==3) || (myOrder[1]==3 && myOrder[2]==0 && myOrder[3]==3) || (myOrder[1]==2 && myOrder[2]==1 && myOrder[3]==2) || (myOrder[1]==5 && myOrder[2]==0 && myOrder[3]==4)|| (myOrder[1]==5 && myOrder[2]==0 && myOrder[3]==3) || length(myOrder[myOrder>4])>0)
        myOrder=c(2,0,2)
    }
    else
      myOrder=c(2,0,2)
  }
  else
  {
    print("WARNING: Could not determine the model !!! (P, d, q) = (2.0.2)!!!")
    myOrder=c(2,0,2)
  }
  
  return (myOrder)
}

calculate_arima_pred=function(timeseries,o,x){
  #datapred=timeseries[(per+n*(day-1)):((ndays*n)+(per-1)+n*(day-1))]
  datapred=timeseries[(1:(x-1))]
  r <- tryCatch((predict(arima(datapred,order=o ), n.ahead =1)$pred), error=function(e) e)
  if (is.numeric(r))
  {if(r<0){r=0}else{r=(r)}}else{r=0}
  return(r)
}


calculate_arima_pred1=function(timeseries,o){
  #datapred=timeseries[(per+n*(day-1)):((ndays*n)+(per-1)+n*(day-1))]
  datapred=timeseries
  r <- tryCatch((predict(arima(datapred,order=o ), n.ahead =1)$pred), error=function(e) e)
  if (is.numeric(r))
  {if(r<0){r=0}else{r=(r)}}else{r=0}
  return(r)
}
pred.arima=function(o,train.ts)
{
  data=rbind(train.ts,test.ts)
  n=nrow(train.ts)
  n.test=nrow(test.ts)
  data.l=lapply(1:n.test, function(x) data[((x)):(n+(x-1)),])
  r=sapply(1:n.test, function(x) calculate_arima_pred1(  data.l[[x]],o))
  return(r)
}

calculate_arima_prediction=function(train.ts, test.ts,H,x){
  
  data=rbind(train.ts,test.ts)
  n=nrow(train.ts)
  n.test=nrow(test.ts)
  
  n.times=n.test/H
  
  ts_train= data[(1+H*(x-1)):(n+H*(x-1)),]
  o=calculaModeloArima(ts_train)
  
  #datapred=data[(1+H*(x-1)):(n+H*(x)),]
  r=sapply(1:H, function(y) calculate_arima_pred(data[(1+H*(x-1)):(n+H*(x-1)+y),],o,(n+y)))
  return(r)
}


calculensemble=function(H,pr,list.pred, real){
  
  rtest=lapply(1: length(list.pred) ,function(i) sapply(1:H, function(y) list.pred[[i]][(pr-H-1+y)]))
  xtest=sapply(1:H, function(y) real[(pr-H-1+y)])
  
  ro=sapply(1: length(list.pred) ,function(i)rmse(rtest[[i]],xtest))
  pos=which(ro==Inf)
  if(length(pos)>0){ro[pos]=rep(1,length(pos))}
  ro.in=1-ro
  gama=sum(ro.in)
  preds=sapply(1: length(list.pred) ,function(i) list.pred[[i]]*ro.in[i])
  r=round( sum(preds[pr,])/gama)
  return(r)
  
}#to  recheck


#pred.ens=sapply(9:672, function(x) calculensemble(H,x,list.pred, real))

#calsmp(672,pred.ens,test.reg$target[9:672])*100
#####

#genrate model lists trained
train.models=function(data.train, formula,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
{
  
  #train models on data train
  
  train.reg= data.train
  
  train.ts=data.frame(data.train.n[,1])
  
  names(train.ts)='ts'
  
  
  #pred.a=unlist(lapply(1:times, function(x)calculate_arima_prediction(train.ts, val.ts,per.arima,x)))
  #svr 
  svr.rbf <-ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[1],kpar ="automatic", prob.model = TRUE)
  #pred.svr.rbf=as.vector(predict(svr.rbf,val.reg, type='response'))
  svr.pol <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[2],kpar ="automatic", prob.model = TRUE)
  #pred.svr.pol=as.vector(predict(svr.pol,val.reg, type='response'))
  
  svr.ln <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  #pred.svr.ln=as.vector(predict(svr.ln,val.reg, type='response'))
  
  
  svr.lap <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[4],kpar ="automatic", prob.model = TRUE)
  #pred.svr.lap=as.vector(predict(svr.lap,val.reg, type='response'))
  
  
  mars.1.7=mars(x=train.reg[,2: ncol(train.reg)], y=train.reg[,1], degree = 1, nk=7)
  
  mars.1.15=mars(x=train.reg[,2: ncol(train.reg)], y=train.reg[,1], degree = 1, nk=15)
  
  
  mars.3.7=mars(x=train.reg[,2: ncol(train.reg)], y=train.reg[,1], degree = 3, nk=7)
  
  mars.3.15=mars(x=train.reg[,2: ncol(train.reg)], y=train.reg[,1], degree = 3, nk=15)
  
  #pred.mars=as.vector(predict(mars,val.reg[,2: ncol(train.reg)]))
  
  
  gbm.10=gbm(formula,data=train.reg,n.trees=10)
  #pred.gbm=as.vector(predict(gbm,val.reg, type='response',n.trees=10))
  
  rf.100 <- randomForest::randomForest(formula,data=train.reg,ntree=100)
  #pred.rf.100=as.vector(predict(rf.100,val.reg, type='response'))
  
  
  rf.250 <- randomForest::randomForest(formula,data=train.reg,ntree=250)
  #pred.rf.250=as.vector(predict(rf.250,val.reg, type='response'))
  
  
  rf.500 <- randomForest(formula,data=train.reg,ntree=500)
  #pred.rf.500=as.vector(predict(rf.500,val.reg, type='response'))
  
  
  
  ppr.2= ppr(formula,data =train.reg, nterms = 2, max.terms = 5)
  #pred.ppr.2=as.vector(predict(ppr.2,val.reg, type='response'))
  
  
  ppr.5= ppr(formula,data =train.reg, nterms =5, max.terms = 5)
  #pred.ppr.5=as.vector(predict(ppr.5,val.reg, type='response'))
  
  
  # mlp.5=mlp(x=train.reg[,2: ncol(train.reg)], y=train.reg[,1],size=c(5),maxit=1000)
  mlp.3= monmlp.fit(x=as.matrix(train.reg[,2: ncol(train.reg)]), y=as.matrix(train.reg[,1]), hidden1=3, n.ensemble=1, bag=F,silent=T)
  #pred.mlp.5=as.vector(predict(mlp.5,as.matrix(test.reg[,2:6])))
  
  #mlp.10=mlp(x=train.reg[,2: ncol(train.reg)], y=train.reg[,1],size=c(10))
  #pred.mlp.10=as.vector(predict(mlp.10,val.reg[,1:5]))
  mlp.5= monmlp.fit(x=as.matrix(train.reg[,2: ncol(train.reg)]), y=as.matrix(train.reg[,1]), hidden1=5, n.ensemble=1, bag=F,silent=T)
  
  
  mlp.7= monmlp.fit(x=as.matrix(train.reg[,2: ncol(train.reg)]), y=as.matrix(train.reg[,1]), hidden1=7, n.ensemble=1, bag=F,silent=T)
  
  #mlp.25=mlp(x=train.reg[,2: ncol(train.reg)], y=train.reg[,1],size=c(25))
  #pred.mlp.25=as.vector(predict(mlp.25,val.reg[,1:5]))
  mlp.10= monmlp.fit(x=as.matrix(train.reg[,2: ncol(train.reg)]), y=as.matrix(train.reg[,1]), hidden1=10, n.ensemble=1, bag=F,silent=T)
  
  mlp.15= monmlp.fit(x=as.matrix(train.reg[,2: ncol(train.reg)]), y=as.matrix(train.reg[,1]), hidden1=15, n.ensemble=1, bag=F,silent=T)
  
  mlp.25= monmlp.fit(x=as.matrix(train.reg[,2: ncol(train.reg)]), y=as.matrix(train.reg[,1]), hidden1=25, n.ensemble=1, bag=F,silent=T)
  
  
  gp.rbf <- gausspr(formula, data =train.reg,type = "regression", kernel = ker[1],kpar ="automatic", prob.model = TRUE)
  #pred.gp.rbf=as.vector(predict(gp.rbf,val.reg, type='response'))
  
  
  gp.pol<- gausspr(formula, data =train.reg,type = "regression", kernel = ker[2],kpar ="automatic", prob.model = TRUE)
  #pred.gp.pol=as.vector(predict(gp.pol,val.reg, type='response'))
  
  
  gp.ln <- gausspr(formula, data =train.reg,type = "regression", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  #pred.gp.ln=as.vector(predict(gp.ln,val.reg, type='response'))
  
  
  gp.lap <- gausspr(formula, data =train.reg,type = "regression", kernel = ker[4],kpar ="automatic", prob.model = TRUE)
  #pred.gp.lap=as.vector(predict(gp.lap,val.reg, type='response'))
  
  
  
  pls.sim = plsr(formula, data=train.reg, validation="CV",method ='simpls')
  
  #pred.pls.sim=as.vector(predict(pls.sim,val.reg, type='response',ncomp=1))
  
  pls.ker=plsr(formula, data=train.reg,method ='widekernelpls')
  #pred.pls.ker=as.vector(predict(pls.ker,val.reg, type='response',ncomp=1))
  
  
  
  pcmr=pcr(formula, data=train.reg)
  #pred.pcmr=as.vector(predict(pcmr,val.reg, type='response',ncomp=1))
  
  
  # dt
  
  
  dt<- rpart(formula, method="anova", data=train.reg )
  
  
  gbm.5=gbm(formula,data=train.reg,n.trees=5)
  
  output.model=list(pls.sim,pls.ker,pcmr,gbm.10,rf.100
                    ,rf.250,rf.500,svr.lap
                    ,svr.pol,svr.rbf,svr.ln,gp.lap,gp.pol
                    ,gp.rbf,gp.ln,mlp.3,mlp.5,mlp.7, mlp.10,mlp.15,mlp.25,ppr.2,ppr.5,
                    mars.1.7,mars.1.15,mars.3.7,mars.3.15,dt,gbm.5)
  
  
  return((output.model))
  
}



predict.models=function(models,data.test,formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
{
  
  val.reg=data.test
  val.ts=data.frame(val.reg[,1])
  names(val.ts)='ts'
  #times=nrow(val.ts)/per.arima
  st=Sys.time()
  
  #pls.sim = plsr(formula, data=train.reg, validation="CV",method ='simpls')
  
  pred.pls.sim=as.vector(predict(models[[1]],val.reg, type='response',ncomp=1))
  
  # pls.ker=plsr(formula, data=train.reg,method ='widekernelpls')
  pred.pls.ker=as.vector(predict(models[[2]],val.reg, type='response',ncomp=1))
  
  
  
  #pcmr=pcr(formula, data=train.reg)
  pred.pcmr=as.vector(predict(models[[3]],val.reg, type='response',ncomp=1))
  #gbm=gbm(formula,data=train.reg,n.trees=10)
  pred.gbm.10=as.vector(predict(models[[4]],val.reg, type='response',n.trees=10))
  
  #rf.100 <- randomForest::randomForest(formula,data=train.reg,ntree=100)
  pred.rf.100=as.vector(predict(models[[5]],val.reg, type='response'))
  
  
  #rf.250 <- randomForest::randomForest(formula,data=train.reg,ntree=250)
  pred.rf.250=as.vector(predict(models[[6]],val.reg, type='response'))
  
  
  #rf.500 <- randomForest(formula,data=train.reg,ntree=500)
  pred.rf.500=as.vector(predict(models[[7]],val.reg, type='response'))
  
  #svr.lap <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[4],kpar ="automatic", prob.model = TRUE)
  pred.svr.lap=as.vector(predict(models[[8]],val.reg, type='response'))
  
  
  pred.svr.pol=as.vector(predict(models[[9]],val.reg, type='response'))
  pred.svr.rbf=as.vector(predict(models[[10]],val.reg, type='response'))
  #svr.pol <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[2],kpar ="automatic", prob.model = TRUE)
  
  
  #svr.ln <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  pred.svr.ln=as.vector(predict(models[[11]],val.reg, type='response'))
  
  
  
  pred.gp.lap=as.vector(predict(models[[12]],val.reg, type='response'))
  
  #gp.pol<- gausspr(formula, data =train.reg,type = "regression", kernel = ker[2],kpar ="automatic", prob.model = TRUE)
  pred.gp.pol=as.vector(predict(models[[13]],val.reg, type='response'))
  
  
  
  
  #gp.rbf <- gausspr(formula, data =train.reg,type = "regression", kernel = ker[1],kpar ="automatic", prob.model = TRUE)
  pred.gp.rbf=as.vector(predict(models[[14]],val.reg, type='response'))
  
  #gp.ln <- gausspr(formula, data =train.reg,type = "regression", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  pred.gp.ln=as.vector(predict(models[[15]],val.reg, type='response'))
  
  
  # gp.lap <- gausspr(formula, data =train.reg,type = "regression", kernel = ker[4],kpar ="automatic", prob.model = TRUE)
  
  
  #mlp.5=mlp(train.reg[,1:5],train.reg[,6],size=c(5))
  pred.mlp.3= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[16]]))
  
  #=as.vector(predict(models[[16]],val.reg[,2: ncol(train.reg)]))
  pred.mlp.5= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[17]]))
  
  
  pred.mlp.7= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[18]]))
  
  
  pred.mlp.10= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[19]]))
  
  #=as.vector(predict(models[[16]],val.reg[,2: ncol(train.reg)]))
  pred.mlp.15= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[20]]))
  
  
  pred.mlp.25= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[21]]))
  
  
  #mlp.10=mlp(train.reg[,1:5],train.reg[,6],size=c(10))
  
  
  #mlp.25=mlp(train.reg[,1:5],train.reg[,6],size=c(25))
  
  
  
  
  #pred.a=unlist(lapply(1:times, function(x)calculate_arima_prediction(train.ts, val.ts,per.arima,x)))
  #ed=Sys.time()
  #svr 
  #svr.rbf <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[1],kpar ="automatic", prob.model = TRUE)
  
  
  
  
  #ppr.2= ppr(formula,data =train.reg, nterms = 2, max.terms = 5)
  pred.ppr.2=as.vector(predict(models[[22]],val.reg, type='response'))
  
  
  #ppr.5= ppr(formula,data =train.reg, nterms =5, max.terms = 5)
  pred.ppr.5=as.vector(predict(models[[23]],val.reg, type='response'))
  
  
  # mars=mars(x=train.reg[,2: ncol(train.reg)], y=train.reg[,1])
  
  pred.mars.1.7=as.vector(predict(models[[24]],val.reg[,2: ncol(train.reg)]))
  
  pred.mars.1.15=as.vector(predict(models[[25]],val.reg[,2: ncol(train.reg)]))
  
  
  pred.mars.3.7=as.vector(predict(models[[26]],val.reg[,2: ncol(train.reg)]))
  
  pred.mars.3.15=as.vector(predict(models[[27]],val.reg[,2: ncol(train.reg)]))
  
  
  
  pred.dt=predict(models[[28]],val.reg)
  
  pred.gbm.5=as.vector(predict(models[[29]],val.reg, type='response',n.trees=5))
  
  
  output.m=data.frame(val.ts,pred.pls.sim,pred.pls.ker,pred.pcmr,pred.gbm.10,pred.rf.100
                      ,pred.rf.250,pred.rf.500,pred.svr.lap
                      ,pred.svr.pol,pred.svr.rbf,pred.svr.ln,pred.gp.lap,pred.gp.pol
                      ,pred.gp.rbf,pred.gp.ln,pred.mlp.3,pred.mlp.7,pred.mlp.5,pred.mlp.10,pred.mlp.15,pred.mlp.25,pred.ppr.2,pred.ppr.5,
                      pred.mars.1.7,pred.mars.1.15,pred.mars.3.7,pred.mars.3.15,pred.dt,pred.gbm.5)
  
  
  ed=Sys.time()
  return((output.m))
  
}

drift.topk.models.input=function(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
{
  
  n=nrow(data.train)
  data.val=data.train[((n-val.length+1):n), ]
  data.train.n=data.train[(1:(n-val.length)), ]
  
  if(t==1){target=data.val[,1]
  #data.train.n=data.train[(1:(n-val.length)), ]
  data.val.n= data.val
  }else if(t>val.length){  
    target=data.test[(t-val.length):(t-1),1]
    #data.train.n=data.train[(t:(n-val.length+t-1)), ]
    data.val.n=data.test[(t-val.length):(t-1),]}else{
      target=c(data.val[((t):val.length),1 ],data.test[1:(t-1),1])
      data.val.n=rbind(data.val[((t):val.length), ],data.test[1:(t-1),])
    }
  
  #train models on data train
  #and compute predictions for the current avialbale set
  
  #output the matrix prediction / target
  
  #arima train and generate prediction on data cal directely
  
  #train.reg= data.train.n
  val.reg= data.val.n
  train.ts=data.frame(data.train.n[,1])
  val.ts=data.frame(data.val.n[,1])
  names(val.ts)='ts'
  #times=nrow(val.ts)/per.arima
  st=Sys.time()
  
  #pls.sim = plsr(formula, data=train.reg, validation="CV",method ='simpls')
  
  #pls.sim = plsr(formula, data=train.reg, validation="CV",method ='simpls')
  
  pred.pls.sim=as.vector(predict(models[[1]],val.reg, type='response',ncomp=1))
  
  # pls.ker=plsr(formula, data=train.reg,method ='widekernelpls')
  pred.pls.ker=as.vector(predict(models[[2]],val.reg, type='response',ncomp=1))
  
  
  
  #pcmr=pcr(formula, data=train.reg)
  pred.pcmr=as.vector(predict(models[[3]],val.reg, type='response',ncomp=1))
  #gbm=gbm(formula,data=train.reg,n.trees=10)
  pred.gbm.10=as.vector(predict(models[[4]],val.reg, type='response',n.trees=10))
  
  #rf.100 <- randomForest::randomForest(formula,data=train.reg,ntree=100)
  pred.rf.100=as.vector(predict(models[[5]],val.reg, type='response'))
  
  
  #rf.250 <- randomForest::randomForest(formula,data=train.reg,ntree=250)
  pred.rf.250=as.vector(predict(models[[6]],val.reg, type='response'))
  
  
  #rf.500 <- randomForest(formula,data=train.reg,ntree=500)
  pred.rf.500=as.vector(predict(models[[7]],val.reg, type='response'))
  
  #svr.lap <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[4],kpar ="automatic", prob.model = TRUE)
  pred.svr.lap=as.vector(predict(models[[8]],val.reg, type='response'))
  
  
  pred.svr.pol=as.vector(predict(models[[9]],val.reg, type='response'))
  pred.svr.rbf=as.vector(predict(models[[10]],val.reg, type='response'))
  #svr.pol <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[2],kpar ="automatic", prob.model = TRUE)
  
  
  #svr.ln <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  pred.svr.ln=as.vector(predict(models[[11]],val.reg, type='response'))
  
  
  
  pred.gp.lap=as.vector(predict(models[[12]],val.reg, type='response'))
  
  #gp.pol<- gausspr(formula, data =train.reg,type = "regression", kernel = ker[2],kpar ="automatic", prob.model = TRUE)
  pred.gp.pol=as.vector(predict(models[[13]],val.reg, type='response'))
  
  
  
  
  #gp.rbf <- gausspr(formula, data =train.reg,type = "regression", kernel = ker[1],kpar ="automatic", prob.model = TRUE)
  pred.gp.rbf=as.vector(predict(models[[14]],val.reg, type='response'))
  
  #gp.ln <- gausspr(formula, data =train.reg,type = "regression", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  pred.gp.ln=as.vector(predict(models[[15]],val.reg, type='response'))
  
  
  # gp.lap <- gausspr(formula, data =train.reg,type = "regression", kernel = ker[4],kpar ="automatic", prob.model = TRUE)
  
  
  #mlp.5=mlp(train.reg[,1:5],train.reg[,6],size=c(5))
  pred.mlp.3= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[16]]))
  
  #=as.vector(predict(models[[16]],val.reg[,2: ncol(train.reg)]))
  pred.mlp.5= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[17]]))
  
  
  pred.mlp.7= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[18]]))
  
  
  pred.mlp.10= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[19]]))
  
  #=as.vector(predict(models[[16]],val.reg[,2: ncol(train.reg)]))
  pred.mlp.15= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[20]]))
  
  
  pred.mlp.25= as.vector(monmlp.predict(as.matrix(val.reg[,2: ncol(train.reg)]),models[[21]]))
  
  
  #mlp.10=mlp(train.reg[,1:5],train.reg[,6],size=c(10))
  
  
  #mlp.25=mlp(train.reg[,1:5],train.reg[,6],size=c(25))
  
  
  
  
  #pred.a=unlist(lapply(1:times, function(x)calculate_arima_prediction(train.ts, val.ts,per.arima,x)))
  #ed=Sys.time()
  #svr 
  #svr.rbf <- ksvm(formula, data =train.reg,type = "eps-svr", kernel = ker[1],kpar ="automatic", prob.model = TRUE)
  
  
  
  
  #ppr.2= ppr(formula,data =train.reg, nterms = 2, max.terms = 5)
  pred.ppr.2=as.vector(predict(models[[22]],val.reg, type='response'))
  
  
  #ppr.5= ppr(formula,data =train.reg, nterms =5, max.terms = 5)
  pred.ppr.5=as.vector(predict(models[[23]],val.reg, type='response'))
  
  
  # mars=mars(x=train.reg[,2: ncol(train.reg)], y=train.reg[,1])
  
  pred.mars.1.7=as.vector(predict(models[[24]],val.reg[,2: ncol(train.reg)]))
  
  pred.mars.1.15=as.vector(predict(models[[25]],val.reg[,2: ncol(train.reg)]))
  
  
  pred.mars.3.7=as.vector(predict(models[[26]],val.reg[,2: ncol(train.reg)]))
  
  pred.mars.3.15=as.vector(predict(models[[27]],val.reg[,2: ncol(train.reg)]))
  
  
  
  pred.dt=predict(models[[28]],val.reg)
  
  pred.gbm.5=as.vector(predict(models[[29]],val.reg, type='response',n.trees=5))
  
  
  output.m=data.frame(val.ts,pred.pls.sim,pred.pls.ker,pred.pcmr,pred.gbm.10,pred.rf.100
                      ,pred.rf.250,pred.rf.500,pred.svr.lap
                      ,pred.svr.pol,pred.svr.rbf,pred.svr.ln,pred.gp.lap,pred.gp.pol
                      ,pred.gp.rbf,pred.gp.ln,pred.mlp.3,pred.mlp.7,pred.mlp.5,pred.mlp.10,pred.mlp.15,pred.mlp.25,pred.ppr.2,pred.ppr.5,
                      pred.mars.1.7,pred.mars.1.15,pred.mars.3.7,pred.mars.3.15,pred.dt,pred.gbm.5)
  
  
  ed=Sys.time()
  return((output.m))
}

rnormc1=function(ts1,ts2)
{   
  mt=data.frame(ts1,ts2)
  c=cor(mt, use="complete.obs", method="pearson")
  rn=sqrt((1-c)/2)
  return(rn)
}

topk.model.sel=function(models,data.train, data.test,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),lim,tp1)
{ 
  t=1
  
  n=NULL
  
  updated_selection=list(0)
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  
  m <- cor(output.m, use="complete",method="pearson")
  #m=m[complete.cases(m),]
  k1=abs(m)
  
  k1[!is.finite(k1)] <- 0
  #k1=sqrt((1-k1)/2)
  
  
  mdl.sel=match(names(sort(k1[1,],decreasing = T)[1:tp1]),names(k1[1,]))
  
  v1=which(k1[1, mdl.sel[-1]]>lim)
  if(length(v1)>0) {
    mdl.sel1=c(1, mdl.sel[v1+1])
    #o1= calculordervar( datatrain1,stdsel) 
  }
  
  min_d=min(sapply(2:tp1, function(x)    rnormc1(output.m[,1],output.m[, mdl.sel[x]])[1,2]))
  # min_d=min(sapply(2:tp1[w], function(x) k1[std,stdsel[x]]))
  
  updated_selection[[1]]= mdl.sel1
  n[1]=0
  t=2
  repeat{
    st=Sys.time()
    output.m1=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
    
    ed=Sys.time()
    m1<- cor( output.m1, use="complete",method="pearson")
    
    k2=abs(m1)
    
    k2[!is.finite(k2)] <- 0
    
    #k2=sqrt((1-k2)/2)
    
    mdl.sel=match(names(sort(k2[1,],decreasing = T)[1:tp1]),names(k2[1,]))
    
    l_cor=min(sapply(2:tp1, function(x) rnormc1(output.m1[,1],output.m1[,mdl.sel[x]])[1,2]))
    
    
    dd=min_d-l_cor
    if(is.na(dd)==T){dd=0}
    
    
    if(abs(dd)>sqrt(log(1/0.95)/(2*nrow(output.m1)))){
      #n=n+1
      mdl.sel1= mdl.sel
      # min_d=min(sapply(2:tp1[w], function(x) k2[std,stdsel[x]]))
      min_d=l_cor
      
      
      
      v=which(k2[1,mdl.sel[-1]]>lim)
      if(length(v)>0) {
        mdl.sel1=c(1, mdl.sel[v+1])
        
      }
      updated_selection[[t]]= mdl.sel1
      n[t]=1
    }else{
      updated_selection[[t]]= mdl.sel1
      n[t]=0
    }
    t=t+1
    if(t>nrow(test.reg)){break}
  }
  return(list("alarm"=n, "models.sel"=updated_selection))
}

ens.step=function(data.test,predictions.table,updated_selection1,t,H)
{  
  if(t<(H+1)){pred=predictions.table[t,updated_selection1[[t]][-1][1]]}else{  
    cols=updated_selection1[[t]][-1]
    
    pred.list=lapply(1:length(cols), function(x)predictions.table[,cols[x]])
    
    pred=calculensemble(H,t,pred.list,predictions.table[,1])
    
  }
  
  return(pred)
}

ens.step.avg=function(data.test,predictions.table,updated_selection1,t,H)
{  
  if(t<(H+1)){pred=predictions.table[t,updated_selection1[[t]][-1][1]]}else{  
    cols=updated_selection1[[t]][-1]
    
    pred.list=lapply(1:length(cols), function(x)predictions.table[,cols[x]])
    
    preds=sapply(1:length( pred.list), function(x) pred.list[[x]])
    pred=mean(sapply(1:length( pred.list), function(x) preds[t,x]))
  }
  
  return(pred)
}



ens.top.pred=function(data.test,predictions.table,updated_selection1,H)
{
  
  pred.ens=sapply(1:nrow(data.test),function(t) ens.step(data.test,predictions.table,updated_selection1,t,H))
  
  
  return(pred.ens)
}

ens.top.pred.avg=function(data.test,predictions.table,updated_selection1,H)
{
  
  pred.ens=sapply(1:nrow(data.test),function(t) ens.step.avg(data.test,predictions.table,updated_selection1,t,H))
  
  
 
  return(pred.ens)
}


ens.step.all=function(data.test,predictions.table,t,H)
{  
  if(t<(H+1)){pred=predictions.table[t,2]}else{  
    pred.list=lapply(2:ncol(predictions.table), function(x)predictions.table[,x])
    
    pred=calculensemble(H,t,pred.list,predictions.table[,1])}
  
  
  
  return(pred)
}

ens.avg.all=function(data.test,predictions.table,t)
{  
  predictions.table1=predictions.table
  pred.list=lapply(2:ncol(predictions.table1), function(x)predictions.table1[,x])
  preds=sapply(1:length( pred.list), function(x) pred.list[[x]])
  pred=mean(sapply(1:length( pred.list), function(x) preds[t,x]))
  
  
  
  return(pred)
}

repr_fea_extract <- function(x) {
  return(c(mean(x), median(x), max(x), min(x), sd(x)))
}


norm_max <- function(x) {
  return(x/max(x))
}


cluster.models=function(sel.models,k.max)
{ 
  # Use verbose = FALSE to hide computing progression.
  e=  tryCatch({
    
    m.cl=as.matrix(t(sel.models))
    opt.k.search=fviz_nbclust(t(sel.models), kmeans, nstart = 25,  method = "silhouette", nboot =50,k.max)+
      labs(subtitle = "silhouette")
    opt.k=which(opt.k.search$data[,2]==max(opt.k.search$data[,2]))
  }, error=function(e){return(e=0)})
  #call methods
  if(e==0){opt.k=5}else{opt.k=e}
  C1.k<- kmeans(t(sel.models),centers = opt.k)
  
  C1.k$cluster
  names.rep=NULL
  
  #computes representative
  for (i in 1:opt.k){
    if(length(which(C1.k$cluster==i))==1)
    { names.rep[i] <- names(which(C1.k$cluster==i)) }else{ 
      rowsum <- rowSums(abs(t(sel.models)[which(C1.k$cluster==i),] - C1.k$centers[i,])) 
      names.rep[i] <-(names(which.min(rowsum)))}
  }
  
  return((names.rep))
}



#keep same models do   not retrain

#add how to generate update selection and n


# recompute clusters once n takes 1 as values

compute.cluster=function(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection,tp1)
{ 
  
  selection=updated_selection$models.sel
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  sel.models.res=abs(output.m[,selection[[t]][-1]]-output.m[,selection[[t]][1]])
  
  sel.models=as.data.frame(output.m[,selection[[t]][-1]])
  names(sel.models)=colnames(output.m[,])[selection[[t]][-1]]
  
  if(ncol(sel.models)<7)
  { names.rep=names(sel.models) }else{ 
    names.rep=cluster.models(sel.models,(ncol(sel.models)-2))}
  #position in the models list
  pos.models=c(1,match(names.rep,names(output.m)))
  
  return(list("names of selected models"=names.rep[[1]],"list.position"=pos.models,  "number of models"=length(pos.models)))
  
}





update.cluster=function(models,data.train, data.test,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection,tp1)
{
  cluster.res=list(0)
  n=updated_selection$alarm
  
  alarm=which(n==1)
  
  cluster.res[[1]]=compute.cluster(models,data.train, data.test,1,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection,tp1)
  for(t in 2:(alarm[1]-1))
  {
    cluster.res[[t]]=cluster.res[[1]]
    
  }
  for(i in 1:(length(alarm)-1))
  {
    cluster.res[[alarm[i]]]=compute.cluster(models,data.train, data.test,alarm[i],val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection,tp1)
    
    for (j in 1:(alarm[(i+1)]-alarm[(i)]-1))  
    {
      cluster.res[[(alarm[i]+j)]]=cluster.res[[alarm[i]]]
      
    }
    
  }
  cluster.res[[alarm[length(alarm)]]]=compute.cluster(models,data.train, data.test,alarm[length(alarm)],val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection,tp1)
  
  for (j in 1:(nrow(data.test)-alarm[length(alarm)]))  
  {
    cluster.res[[(alarm[length(alarm)]+j)]]=cluster.res[[alarm[i]]]
    
  }
  
  
  return(cluster.res)
}


update.cluster1=function(models,data.train, data.test,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection,tp1)
{
  n=updated_selection$alarm
  
  alarm=which(n==1)
  
  cluster.res1=lapply(1:nrow(data.test), function(t) compute.cluster(models,data.train, data.test,t,val.length, formula, per.arima
                                                                     ,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection,tp1))
  
  return(cluster.res1)
}


compute.cluster.tp=function(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection)
{ 
  
  selection=updated_selection$models.sel
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  sel.models.res=abs(output.m[,selection[[t]][-1]]-output.m[,selection[[t]][1]])
  
  sel.models=as.data.frame(output.m[,selection[[t]][-1]])
  names(sel.models)=colnames(output.m[,])[selection[[t]][-1]]
  
  k.max=round(ncol(sel.models)/2)-1
  
  
  if(ncol(sel.models)<4){names.rep=names(sel.models)}else{ 
    
    
    a <- tryCatch({a=otrimle(data= t(sel.models), G=k.max, ncores=1)}, error = function(e) 0)
    
    if(is.numeric(a)  ){  
      G=k.max
      while(is.numeric(a) ){ 
        a <- tryCatch({
          a <- otrimle(data= t(sel.models), G=k.max, ncores=1,initial = sample(0:(G),  size=ncol(sel.models), replace=TRUE))
        }, error = function(e) 0)
      }
    }
    cl=a$cluster
    
    names(cl)=names(sel.models)
    
    if(k.max==1){  names.rep=names(which(cl==1))}else{ 
      vec=c(1:k.max)
      d=sapply(1:k.max, function(x) which(cl==vec[x])[1])
      d=d[!is.na(d)]
      names.rep=names(d)
      
      
    }}
  
  pos.models=c(match(names.rep,names(output.m))) 
  
  pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
  
  pred=calculensemble(H,t,pred.list,predictions.table[,1])
  
  return(pred)
}

compute.cluster.tp.all=function(models,data.train,data.test,val.length,H, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection)
{ 
  
  pred=NULL
  tr.model= list(0)
  pred[1:(H)]=predictions.table[1:(H),1]
  for(t in 1:H){
    alarm=updated_selection$alarm
    
    selection=updated_selection$models.sel
    
    output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
    sel.models.res=abs(output.m[,selection[[t]][-1]]-output.m[,selection[[t]][1]])
    
    sel.models=as.data.frame(output.m[,selection[[t]][-1]])
    names(sel.models)=colnames(output.m[,])[selection[[t]][-1]]
    
    k.max=round(ncol(sel.models)/2)-1
    
    
    if(ncol(sel.models)<4){names.rep=names(sel.models)}else{ 
      
      
      a <- tryCatch({a=otrimle(data= t(sel.models), G=k.max, ncores=1)}, error = function(e) 0)
      
      if(is.numeric(a)  ){  
        G=k.max
        while(is.numeric(a) ){ 
          a <- tryCatch({
            a <- otrimle(data= t(sel.models), G=k.max, ncores=1,initial = sample(0:(G),  size=ncol(sel.models), replace=TRUE))
          }, error = function(e) 0)
        }
      }
      cl=a$cluster
      
      names(cl)=names(sel.models)
      
      if(k.max==1){  names.rep=names(which(cl==1))}else{ 
        vec=c(1:k.max)
        d=sapply(1:k.max, function(x) which(cl==vec[x])[1])
        d=d[!is.na(d)]
        names.rep=names(d)
        
        
      }}
    
    pos.models=c(match(names.rep,names(output.m))) 
    
    tr.model[[t]]=pos.models
    
  }
  t=(H+1)
  
  
  alarm=updated_selection$alarm
  
  selection=updated_selection$models.sel
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  sel.models.res=abs(output.m[,selection[[t]][-1]]-output.m[,selection[[t]][1]])
  
  sel.models=as.data.frame(output.m[,selection[[t]][-1]])
  names(sel.models)=colnames(output.m[,])[selection[[t]][-1]]
  
  k.max=round(ncol(sel.models)/2)-1
  
  
  if(ncol(sel.models)<4){names.rep=names(sel.models)}else{ 
    
    
    a <- tryCatch({a=otrimle(data= t(sel.models), G=k.max, ncores=1)}, error = function(e) 0)
    
    if(is.numeric(a)  ){  
      G=k.max
      while(is.numeric(a) ){ 
        a <- tryCatch({
          a <- otrimle(data= t(sel.models), G=k.max, ncores=1,initial = sample(0:(G),  size=ncol(sel.models), replace=TRUE))
        }, error = function(e) 0)
      }
    }
    cl=a$cluster
    
    names(cl)=names(sel.models)
    
    if(k.max==1){  names.rep=names(which(cl==1))}else{ 
      vec=c(1:k.max)
      d=sapply(1:k.max, function(x) which(cl==vec[x])[1])
      d=d[!is.na(d)]
      names.rep=names(d)
      
      
    }}
  
  pos.models=c(match(names.rep,names(output.m))) 
  
  tr.model[[t]]=pos.models
  
  pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
  
  pred[t]=calculensemble(H,t,pred.list,predictions.table[,1])
  t=t+1
  repeat{
    if(alarm[t]==1) {
      output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
      sel.models.res=abs(output.m[,selection[[t]][-1]]-output.m[,selection[[t]][1]])
      
      sel.models=as.data.frame(output.m[,selection[[t]][-1]])
      names(sel.models)=colnames(output.m[,])[selection[[t]][-1]]
      
      k.max=round(ncol(sel.models)/2)-1
      
      
      if(ncol(sel.models)<4){names.rep=names(sel.models)}else{ 
        
        
        a <- tryCatch({a=otrimle(data= t(sel.models), G=k.max, ncores=1)}, error = function(e) 0)
        
        if(is.numeric(a)  ){  
          G=k.max
          while(is.numeric(a) ){ 
            a <- tryCatch({
              a <- otrimle(data= t(sel.models), G=k.max, ncores=1,initial = sample(0:(G),  size=ncol(sel.models), replace=TRUE))
            }, error = function(e) 0)
          }
        }
        cl=a$cluster
        
        names(cl)=names(sel.models)
        
        if(k.max==1){  names.rep=names(which(cl==1))}else{ 
          vec=c(1:k.max)
          d=sapply(1:k.max, function(x) which(cl==vec[x])[1])
          d=d[!is.na(d)]
          names.rep=names(d)
          
          
        }}
      
      pos.models=c(match(names.rep,names(output.m))) 
      
      tr.model[[t]]=pos.models
      
      pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
      
      pred[t]=calculensemble(H,t,pred.list,predictions.table[,1])
      
      
      
    }else{
      
      pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
      
      pred[t]=calculensemble(H,t,pred.list,predictions.table[,1]) 
      
      tr.model[[t]]=pos.models
      
    }
    
    
    t=t+1
    
    if(t>nrow(data.test)){break}
    
  }
  
  
  plot(data.test$target, type='l')
  lines(pred, col='red')
  
  print( rmse(data.test$target,pred))
  return(list(pred,tr.model))
}


compute.cluster.or=function(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
{ 
  
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  
  sel.models=output.m[,-c(1,17:19)]
  
  
  k.max=round(ncol(sel.models)/2)-1
  
  a <- tryCatch({a <- otrimle(data= t(sel.models), G=k.max, ncores=1,initial = sample(0:(k.max),  size=ncol(sel.models), replace=TRUE))}, error = function(e) 0)
  
  
  
  while(is.numeric(a) ){ 
    a <- tryCatch({a <- otrimle(data= t(sel.models), G=k.max, ncores=1,initial = sample(0:(k.max),  size=ncol(sel.models), replace=TRUE))}, error = function(e) 0)
  }
  cl=a$cluster
  
  names(cl)=names(sel.models)
  
  if(k.max==1){  names.rep=names(which(cl==1))}else{ 
    vec=c(1:k.max)
    d=sapply(1:k.max, function(x) which(cl==vec[x])[1])
    d=d[!is.na(d)]
    names.rep=names(d)
    
    
  }
  
  pos.models=c(match(names.rep,names(output.m))) 
  
  pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
  
  pred=calculensemble(H,t,pred.list,predictions.table[,1])
  
  return(pred)
}


#pred.cl.or=sapply((H+1):nrow(data.test), function(t)compute.cluster.or(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot")))

#pred.cl.or=c(predictions.table[1:8,1],pred.cl.or)








compute.cluster.or.all=function(models,data.train, data.test,val.length,H, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
{ 
  
  pred=NULL
  pred[1:(H)]=predictions.table[1:(H),1]
  
  t=(H+1)
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  
  sel.models=output.m[,-c(1,17:19)]
  
  
  k.max=round(ncol(sel.models)/2)-1
  
  a <- tryCatch({a <- otrimle(data= t(sel.models), G=k.max, ncores=1,initial = sample(0:(k.max),  size=ncol(sel.models), replace=TRUE))}, error = function(e) 0)
  
  
  
  while(is.numeric(a) ){ 
    a <- tryCatch({a <- otrimle(data= t(sel.models), G=k.max, ncores=1,initial = sample(0:(k.max),  size=ncol(sel.models), replace=TRUE))}, error = function(e) 0)
  }
  cl=a$cluster
  
  names(cl)=names(sel.models)
  
  if(k.max==1){  names.rep=names(which(cl==1))}else{ 
    vec=c(1:k.max)
    d=sapply(1:k.max, function(x) which(cl==vec[x])[1])
    d=d[!is.na(d)]
    names.rep=names(d)
    
    
  }
  
  pos.models=c(match(names.rep,names(output.m))) 
  
  pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
  
  
  repeat{  
    
    pred[t]=calculensemble(H,t,pred.list,predictions.table[,1])
    
    t=t+1
    
    if(t>nrow(data.test)){break}
  }
  
  plot(data.test$target, type='l')
  lines(pred, col='red')
  
  print( rmse(data.test$target,pred))
  
  
  return(pred)
}


library("dtwclust")


compute.cluster.tw.tp=function(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection)
{ 
  
  selection=updated_selection$models.sel
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  sel.models.res=abs(output.m[,selection[[t]][-1]]-output.m[,selection[[t]][1]])
  
  sel.models=as.data.frame(output.m[,selection[[t]][-1]])
  names(sel.models)=colnames(output.m[,])[selection[[t]][-1]]
  
  k.max=round(ncol(sel.models)/2)-1
  
  if( k.max<2){names.rep=names(sel.models)}else{ 
    
    a=tsclust(series =t(sel.models), type = "partitional", k = k.max,
              distance = "dtw_basic")
    
    id=a@control$distmat$id_cent
    cl=a@cluster
    names(cl)=names(sel.models)
    names.rep=names(cl)[id]}
  pos.models=c(match(names.rep,names(output.m))) 
  
  pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
  
  pred=calculensemble(H,t,pred.list,predictions.table[,1])
  
  return(pred)
}



compute.cluster.tw.tp.all=function(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection)
{ 
  pred=NULL
  pred[1:(H)]=predictions.table[1:(H),1]
  
  t=(H+1)
  
  
  alarm=updated_selection$alarm
  
  
  
  selection=updated_selection$models.sel
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  sel.models.res=abs(output.m[,selection[[t]][-1]]-output.m[,selection[[t]][1]])
  
  sel.models=as.data.frame(output.m[,selection[[t]][-1]])
  names(sel.models)=colnames(output.m[,])[selection[[t]][-1]]
  
  k.max=round(ncol(sel.models)/2)-1
  
  if( k.max<2){names.rep=names(sel.models)}else{ 
    
    a=tsclust(series =t(sel.models), type = "partitional", k = k.max,
              distance = "dtw_basic")
    
    id=a@control$distmat$id_cent
    cl=a@cluster
    names(cl)=names(sel.models)
    names.rep=names(cl)[id]}
  pos.models=c(match(names.rep,names(output.m))) 
  
  pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
  
  pred[t]=calculensemble(H,t,pred.list,predictions.table[,1])
  
  t=t+1
  repeat{
    if(alarm[t]==1) { output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
    sel.models.res=abs(output.m[,selection[[t]][-1]]-output.m[,selection[[t]][1]])
    
    sel.models=as.data.frame(output.m[,selection[[t]][-1]])
    names(sel.models)=colnames(output.m[,])[selection[[t]][-1]]
    
    k.max=round(ncol(sel.models)/2)-1
    
    if( k.max<2){names.rep=names(sel.models)}else{ 
      
      a=tsclust(series =t(sel.models), type = "partitional", k = k.max,
                distance = "dtw_basic")
      
      id=a@control$distmat$id_cent
      cl=a@cluster
      names(cl)=names(sel.models)
      names.rep=names(cl)[id]}
    pos.models=c(match(names.rep,names(output.m))) 
    
    pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
    
    pred[t]=calculensemble(H,t,pred.list,predictions.table[,1])}else{
      
      pred.list=lapply(1:length(pos.models), function(x)predictions.table[,pos.models[x]])
      
      pred[t]=calculensemble(H,t,pred.list,predictions.table[,1])
      
    }
    
    t=t+1
    
    if(t>nrow(data.test)){break}
    
  }
  
  return(pred)
}


library(Cubist)


library(glmnet)

compute.cl.tp.stac=function(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),pos.models,updated_selection)
{ 
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  
  
  sel.models=(output.m[,c(1,pos.models[[t]])])
  
  x= sel.models[,2:ncol(sel.models)]
  y= sel.models[,1]
  train.st=sel.models
  formula1=ts~.
  
  test.st=predictions.table[t,c(1,pos.models[[t]])]
  names(test.st)=names(train.st)
  x.test= test.st[,2:ncol(sel.models)]
  y.test= test.st[,1]
  
  
  
  rf.st <- randomForest::randomForest(ts~.,data=train.st,ntree=10)
  
  pred.rfst=predict(rf.st, test.st, type="response")
  
  
  rbr.st=cubist(x, y, committees = 1, control = cubistControl(), weights = NULL)
  pred.rbrst=predict(rbr.st, x.test, type="response")
  
  
  
  
  glm.st=glmnet(as.matrix(x),y, alpha = 0.2, nlambda = 20)
  
  pred.glmst=predict(glm.st, newx = as.matrix(x.test), type = "response")[10]
  
  
  
  svr.st <- ksvm(formula1, data =train.st,type = "eps-svr", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  pred.svrst=predict(svr.st, test.st, type="response")
  
  
  mars.st=mars(x=train.st[,2: ncol(train.st)], y=train.st[,1], degree =1, nk=7)
  pred.marsst=predict(mars.st,x.test, type="response")
  
  
  gbm.st=gbm(formula1,data=train.st,n.trees=10)
  
  pred.gbmst=predict(gbm.st,x.test, type="response",n.trees=10)
  
  
  ppr.st= ppr(formula1,data =train.st, nterms =5, max.terms = 5)
  
  pred.pprst=predict(ppr.st,x.test, type="response")
  
  
  mlp.st= monmlp.fit(x=as.matrix(train.st[,2: ncol(train.st)]), y=as.matrix(train.st[,1]), hidden1=5, n.ensemble=1, bag=F,silent=T)
  pred.mlpst= as.vector(monmlp.predict(as.matrix(test.st[,2: ncol(test.st)]),mlp.st))
  
  
  gp.st <- gausspr(formula1, data =train.st,type = "regression", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  
  pred.gpst=as.vector(predict(gp.st,test.st, type='response'))
  
  
  
  
  pls.st = plsr(formula1, data=train.st, validation="CV",method ='simpls')
  
  pred.plsst=as.vector(predict(pls.st,test.st, type='response',ncomp=1))
  
  
  dt.st<- rpart(formula1, method="anova", data=train.st )
  
  pred.dtst=predict(dt.st,test.st)
  
  vec.pred=c(pred.rfst,pred.rbrst, pred.glmst,pred.svrst,pred.marsst
             ,pred.gbmst,pred.pprst,pred.mlpst,pred.gpst,pred.plsst,pred.dtst)
  
  mean.pred=mean(vec.pred) 
  
  
  return(c(vec.pred,mean.pred))
  
}



compute.cl.tp.stac.all=function(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),pos.models,updated_selection)
{ 
  t=1
  alarm=updated_selection$alarm
  pred.rfst=NULL
  pred.rbrst=NULL
  pred.glmst=NULL
  pred.svrst=NULL
  pred.marsst=NULL
  
  pred.gbmst=NULL
  pred.pprst=NULL
  pred.mlpst=NULL
  pred.gpst=NULL
  pred.plsst=NULL
  pred.dtst=NULL
  
  mean.pred=NULL
  
  
  
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  
  
  sel.models=(output.m[,c(1,pos.models[[t]])])
  
  x= sel.models[,2:ncol(sel.models)]
  y= sel.models[,1]
  train.st=sel.models
  formula1=ts~.
  
  test.st=predictions.table[t,c(1,pos.models[[t]])]
  names(test.st)=names(train.st)
  x.test= test.st[,2:ncol(sel.models)]
  y.test= test.st[,1]
  
  pos.models.k= pos.models[[t]]
  
  rf.st <- randomForest::randomForest(ts~.,data=train.st,ntree=10)
  
  pred.rfst[t]=predict(rf.st, test.st, type="response")
  
  
  rbr.st=cubist(x, y, committees = 1, control = cubistControl(), weights = NULL)
  pred.rbrst[t]=predict(rbr.st, x.test, type="response")
  
  
  
  
  glm.st=glmnet(as.matrix(x),y, alpha = 0.2, nlambda = 20)
  
  pred.glmst[t]=predict(glm.st, newx = as.matrix(x.test), type = "response")[10]
  
  
  
  svr.st <- ksvm(formula1, data =train.st,type = "eps-svr", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  pred.svrst[t]=predict(svr.st, test.st, type="response")
  
  
  mars.st=mars(x=train.st[,2: ncol(train.st)], y=train.st[,1], degree =1, nk=7)
  pred.marsst[t]=predict(mars.st,x.test, type="response")
  
  
  gbm.st=gbm(formula1,data=train.st,n.trees=2)
  
  pred.gbmst[t]=predict(gbm.st,x.test, type="response",n.trees=10)
  
  
  ppr.st= ppr(formula1,data =train.st, nterms =5, max.terms = 5)
  
  pred.pprst[t]=predict(ppr.st,x.test, type="response")
  
  
  mlp.st= monmlp.fit(x=as.matrix(train.st[,2: ncol(train.st)]), y=as.matrix(train.st[,1]), hidden1=5, n.ensemble=1, bag=F,silent=T)
  pred.mlpst[t]= as.vector(monmlp.predict(as.matrix(test.st[,2: ncol(test.st)]),mlp.st))
  
  
  gp.st <- gausspr(formula1, data =train.st,type = "regression", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
  
  pred.gpst[t]=as.vector(predict(gp.st,test.st, type='response'))
  
  
  
  
  pls.st = plsr(formula1, data=train.st, validation="CV",method ='simpls')
  
  pred.plsst[t]=as.vector(predict(pls.st,test.st, type='response',ncomp=1))
  
  
  dt.st<- rpart(formula1, method="anova", data=train.st )
  
  pred.dtst[t]=predict(dt.st,test.st)
  
  
  vec.pred=c(pred.rfst[t],pred.rbrst[t], pred.glmst[t],pred.svrst[t]
             ,pred.gpst[t],pred.plsst[t])
  
  
  mean.pred[t]=mean(vec.pred) 
  
  t=t+1
  
  repeat{
    
    if(alarm[t]==1){ output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
    
    
    sel.models=(output.m[,c(1,pos.models[[t]])])
    
    x= sel.models[,2:ncol(sel.models)]
    y= sel.models[,1]
    train.st=sel.models
    formula1=ts~.
    
    test.st=predictions.table[t,c(1,pos.models[[t]])]
    names(test.st)=names(train.st)
    x.test= test.st[,2:ncol(sel.models)]
    y.test= test.st[,1]
    
   pos.models.k= pos.models[[t]]
    
    rf.st <- randomForest::randomForest(ts~.,data=train.st,ntree=10)
    
    pred.rfst[t]=predict(rf.st, test.st, type="response")
    
    
    rbr.st=cubist(x, y, committees = 1, control = cubistControl(), weights = NULL)
    pred.rbrst[t]=predict(rbr.st, x.test, type="response")
    
    
    
    
    glm.st=glmnet(as.matrix(x),y, alpha = 0.2, nlambda = 20)
    
    pred.glmst[t]=predict(glm.st, newx = as.matrix(x.test), type = "response")[10]
    
    
    
    svr.st <- ksvm(formula1, data =train.st,type = "eps-svr", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
    pred.svrst[t]=predict(svr.st, test.st, type="response")
    
    
    mars.st=mars(x=train.st[,2: ncol(train.st)], y=train.st[,1], degree =1, nk=7)
    pred.marsst[t]=predict(mars.st,x.test, type="response")
    
    
    gbm.st=gbm(formula1,data=train.st,n.trees=5)
    
    pred.gbmst[t]=predict(gbm.st,x.test, type="response",n.trees=10)
    
    
    ppr.st= ppr(formula1,data =train.st, nterms =5, max.terms = 5)
    
    pred.pprst[t]=predict(ppr.st,x.test, type="response")
    
    
    mlp.st= monmlp.fit(x=as.matrix(train.st[,2: ncol(train.st)]), y=as.matrix(train.st[,1]), hidden1=5, n.ensemble=1, bag=F,silent=T)
    pred.mlpst[t]= as.vector(monmlp.predict(as.matrix(test.st[,2: ncol(test.st)]),mlp.st))
    
    
    gp.st <- gausspr(formula1, data =train.st,type = "regression", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
    
    pred.gpst[t]=as.vector(predict(gp.st,test.st, type='response'))
    
    
    
    
    pls.st = plsr(formula1, data=train.st, validation="CV",method ='simpls')
    
    pred.plsst[t]=as.vector(predict(pls.st,test.st, type='response',ncomp=1))
    
    
    dt.st<- rpart(formula1, method="anova", data=train.st )
    
    pred.dtst[t]=predict(dt.st,test.st)
    
    vec.pred=c(pred.rfst[t],pred.rbrst[t], pred.glmst[t],pred.svrst[t]
               ,pred.gpst[t],pred.plsst[t])
    
    mean.pred[t]=mean(vec.pred) }else{
      
      output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
      
      
      sel.models=(output.m[,c(1,pos.models.k)])
      
      x= sel.models[,2:ncol(sel.models)]
      y= sel.models[,1]
      train.st=sel.models
      formula1=ts~.
      
      test.st=predictions.table[t,c(1,pos.models.k)]
      names(test.st)=names(train.st)
      x.test= test.st[,2:ncol(sel.models)]
      y.test= test.st[,1]
      
      
      
      #rf.st <- randomForest::randomForest(ts~.,data=train.st,ntree=10)
      
      pred.rfst[t]=predict(rf.st, test.st, type="response")
      
      
      #rbr.st=cubist(x, y, committees = 1, control = cubistControl(), weights = NULL)
      pred.rbrst[t]=predict(rbr.st, x.test, type="response")
      
      
      
      
      #glm.st=glmnet(as.matrix(x),y, alpha = 0.2, nlambda = 20)
      
      pred.glmst[t]=predict(glm.st, newx = as.matrix(x.test), type = "response")[10]
      
      
      
      #svr.st <- ksvm(formula1, data =train.st,type = "eps-svr", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
      pred.svrst[t]=predict(svr.st, test.st, type="response")
      
      
      #mars.st=mars(x=train.st[,2: ncol(train.st)], y=train.st[,1], degree =1, nk=7)
      pred.marsst[t]=predict(mars.st,x.test, type="response")
      
      
      #gbm.st=gbm(formula1,data=train.st,n.trees=10)
      
      pred.gbmst[t]=predict(gbm.st,x.test, type="response",n.trees=10)
      
      
      #ppr.st= ppr(formula1,data =train.st, nterms =5, max.terms = 5)
      
      pred.pprst[t]=predict(ppr.st,x.test, type="response")
      
      
      #mlp.st= monmlp.fit(x=as.matrix(train.st[,2: ncol(train.st)]), y=as.matrix(train.st[,1]), hidden1=5, n.ensemble=1, bag=F,silent=T)
      pred.mlpst[t]= as.vector(monmlp.predict(as.matrix(test.st[,2: ncol(test.st)]),mlp.st))
      
      
      #gp.st <- gausspr(formula1, data =train.st,type = "regression", kernel = ker[3],kpar ="automatic", prob.model = TRUE)
      
      pred.gpst[t]=as.vector(predict(gp.st,test.st, type='response'))
      
      
      
      
      #pls.st = plsr(formula1, data=train.st, validation="CV",method ='simpls')
      
      pred.plsst[t]=as.vector(predict(pls.st,test.st, type='response',ncomp=1))
      
      
      #dt.st<- rpart(formula1, method="anova", data=train.st )
      
      pred.dtst[t]=predict(dt.st,test.st)
      
      vec.pred=c(pred.rfst[t],pred.rbrst[t], pred.glmst[t],pred.svrst[t]
                 ,pred.gpst[t],pred.plsst[t])
      
      mean.pred[t]=mean(vec.pred) 
    }
    t=t+1
    
    if(t>nrow(data.test)){break}
  }
  
  
  return(list(pred.rfst,pred.rbrst, pred.glmst,pred.svrst,pred.marsst,pred.gbmst
              ,pred.pprst,pred.mlpst,pred.gpst,pred.plsst,pred.dtst
              , mean.pred ))
  
}
compute.stacking=function(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),pos.models,updated_selection)
{ 
  
  output.m=drift.topk.models.input(models,data.train, data.test,t,val.length, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))
  
  
  sel.models=(output.m[,c(1,pos.models[[t]])])
  
  x= sel.models[,2:ncol(sel.models)]
  y= sel.models[,1]
  train.st=sel.models
  formula1=ts~.
  
  test.st=predictions.table[t,c(1,pos.models[[t]])]
  names(test.st)=names(train.st)
  x.test= test.st[,2:ncol(sel.models)]
  y.test= test.st[,1]
  
  fitControl <- trainControl(
    method = 'cv',                   # k-fold cross validation
    number = 5,                      # number of folds
    savePredictions = 'final'      # saves predictions for optimal tuning parameter
    
  ) 
  
  mlboost=train(formula1,train.st,method='gamboost',  trControl = fitControl)
  
  pred.bst=predict(mlboost,test.st)
  
  
  
  
  
  return(pred.bst)
  
}



bias.variance.dcmp=function(y.hat,y.test)
{
  
  bias = var(((y.hat)-y.test))^2
  variance = mean(var(y.hat))
  
  
  return(c(variance=variance,bias_squared= bias))
}


Wilcoxon.test=function(real,pred.1,pred.2,w)
{
  
  n=length(real)/w
  pred_1=lapply(1:w, function(x)pred.1[(1+(x-1)*n):(n*x)])
  pred_2=lapply(1:w, function(x)pred.2[(1+(x-1)*n):(n*x)])
  real_d=lapply(1:w, function(x)real[(1+(x-1)*n):(n*x)])
  
  #pred_1=pred.1
  #pred_2=pred.2
  #real_d=real
  
  
  #res1=abs(real_d-pred_1)^4
  
  res1=sapply(1:w, function(x)rmse(pred_1[[x]],real_d[[x]]))
  
  res2=sapply(1:w, function(x)rmse(pred_2[[x]],real_d[[x]]))
  
  
  #res2=abs(real_d-pred_2)^4
  
  plot(res1,type='l')
  
  lines(res2,col='red')
  method=c(rep('method1',length(pred_1)),rep('method2',length(pred_2)))
  
  res=c(res1, res2)
  
  data=data.frame(method,res)
  
  res=wilcox.test(res ~ method, data = data, paired = T, alternative = "less")
  
  
  res1=t.test(res ~ method, data =data, alternative = "less") 
  
  p.value=res$p.value
  
  p.value1=res1$p.value
  
  
  return(c(p.value,p.value1))
}


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


path='/Users/saadalla/Desktop/ecml-exp2'
setwd(path)





#data=read.delim("/Users/saadalla/Desktop/ecml-exp2/meteo-data-2.txt",sep=',')


#data=read.csv('/Users/saadalla/Desktop/ecml-exp2/Bike-Sharing-Dataset/hour.csv')
data=read.csv('/Users/saadalla/Documents/taxi-data-sets/porto_new/Demand1.csv')


data=data[9001:11693,]
plot(data[,11],type='l')

data=data.frame(data[1:2693,11])
dataset=data_lag_prep(data,1,5)



train.ts=data.frame(dataset[1:2000,1])
test.ts=data.frame(dataset[2001:2688,1])

names(train.ts)=names(test.ts)="target"

train.reg=dataset[1:2000,]
test.reg=dataset[2001:2688,]

#evaluation

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





target=data.test$target

pred.rf=predictions.table$pred.rf.500

pred.gbm=predictions.table$pred.gbm.10

final.preds=data.frame(target, pred.rf, pred.gbm,pred.ens.all,pred.ens.all.sw,pred.a,pred.ade@y_hat,pred.dets@y_hat,
                       pred.ewa,pred.fs,pred.ogd,pred.mlpol,pred.stacking
                       ,pred.cl.or,pred.top,pred.cl,pred.tw.tp,pred.cl.st, pred.tp.cl.all)


rmse.table=sapply(1:ncol(final.preds), function(x) rmse(final.preds[,1],final.preds[,x] ))
names(rmse.table)=names(final.preds)


mse.table=sapply(1:ncol(final.preds), function(x) mse(final.preds[,1],final.preds[,x] ))
names(mse.table)=names(final.preds)


wilcoxon.table=sapply(1:ncol(final.preds), function(x) Wilcoxon.test(final.preds[,1],final.preds[,18],final.preds[,x],688))


tests=round(t(wilcoxon.table),digits = 4)
rownames(tests)=names(final.preds)



alpha1=sapply(1:ncol(final.preds), function(x) if(tests[x,1]<0.5){1}else{0})
names(alpha1)=names(wilcoxon.table)

alpha2=sapply(1:ncol(final.preds), function(x) if(tests[x,2]<0.5){1}else{0})
names(alpha2)=names(wilcoxon.table)



beta1=sapply(1:ncol(final.preds), function(x) if(tests[x,1]<0.06){1}else{0})
names(beta1)=names(wilcoxon.table)

beta2=sapply(1:ncol(final.preds), function(x) if(tests[x,2]<0.05){1}else{0})
names(beta2)=names(wilcoxon.table)



