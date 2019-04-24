


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