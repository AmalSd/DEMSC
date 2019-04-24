data_lag_prep=function(data,col.indx,k1){  
  # Pre-processig and preparibg data
  reg_data=(t(sapply((k1+1):(dim(data)[1]),function(x) return(data[c(x,(x-k1):(x-1)),col.indx]))))
  colnames(reg_data)= c('target',paste("lag",k1:1,sep=''))
  reg_data=as.data.frame(reg_data)
  return(reg_data)
}
