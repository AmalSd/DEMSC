

# drift aware clustering computation




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


