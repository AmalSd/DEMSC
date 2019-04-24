
#different ensemble framework calculations

# weighted averaged ensemble/ simple average 

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
