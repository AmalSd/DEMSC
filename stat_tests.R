
#evaluation 


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
