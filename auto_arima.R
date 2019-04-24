

#Arima order automatic setting


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


#calculate prediction for one step ahead
calculate_arima_pred=function(timeseries,o,x){
  #datapred=timeseries[(per+n*(day-1)):((ndays*n)+(per-1)+n*(day-1))]
  datapred=timeseries[(1:(x-1))]
  r <- tryCatch((predict(arima(datapred,order=o ), n.ahead =1)$pred), error=function(e) e)
  if (is.numeric(r))
  {if(r<0){r=0}else{r=(r)}}else{r=0}
  return(r)
}

#calculate prediction for the whole test set

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

#calculate prediction for each step with periodic update


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
