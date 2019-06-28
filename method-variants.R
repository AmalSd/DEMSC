


H=10# size of the sliding window for computing ensembles candidates error (see the paper)

# sliding window ensemble of the top base models
pred.top=ens.top.pred(data.test,predictions.table,updated_selection1,H)



# sliding window ensemble of all the models in the pool

pred.ens.all.sw=sapply(1:nrow(data.test),function(t)ens.step.all(data.test,predictions.table,t,H))

# average ensemble of all the models in the pool

pred.ens.all=sapply(1:nrow(data.test),function(t)ens.avg.all(data.test,predictions.table,t))



rmse(data.test$target,pred.top)
rmse(data.test$target,pred.ens.all)

rmse(data.test$target,pred.ens.all.sw)



# top model clustering using K-means and then aggregation in a sliding window ensemble
#Clusters are updated with each drift detection
val.length1=10
cluster.res=update.cluster.kmeans(models,data.train, data.test,val.length1, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection,tp1)

pos.cl=lapply(1:length(cluster.res) ,function(x) if(length(cluster.res[[x]]$list.position)>=4){
  cluster.res[[x]]$list.position[1:4]
}else{c(cluster.res[[x]]$list.position[1],rep(cluster.res[[x]]$list.position[2],3))})


pred.cl.kmeans=ens.top.pred(data.test,predictions.table,pos.cl,H)


rmse(data.test$target,pred.cl.kmeans)




# original model pool clustering using IMLEC and then aggregation in a sliding window ensemble 

val.length1=10

pred.cl.or=compute.cluster.or.all(models,data.train, data.test,val.length1,H, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"))



rmse(data.test$target,pred.cl.or)



# top model clustering using dynamic time wrapping  and then aggregation in a sliding window ensemble 
# Clusters are updated with each drift detection

val.length1=10
pred.tw.tp=compute.cluster.tw.tp.all(models,data.train, data.test,t,val.length1, formula, per.arima,ker=c("rbfdot" ,"polydot" ,"vanilladot", "laplacedot"),updated_selection)

rmse(data.test$target,pred.tw.tp)



