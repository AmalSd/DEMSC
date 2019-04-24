

#train models and generate base-line predictions


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
