



rnormc1=function(ts1,ts2)
{   
  mt=data.frame(ts1,ts2)
  c=cor(mt, use="complete.obs", method="pearson")
  rn=sqrt((1-c)/2)
  return(rn)
}


repr_fea_extract <- function(x) {
  return(c(mean(x), median(x), max(x), min(x), sd(x)))
}


norm_max <- function(x) {
  return(x/max(x))
}

