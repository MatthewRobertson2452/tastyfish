

#' Plot the functional response curve and associated data
#'
#' @param rep The report object from the TMB model
#' @param data The datalist used for the TMB model
#'
#' @return plot
#' @export
#'
#' @examples plot_curve(rep=rep, data=data)
plot_curve<-function(rep, data){
  new_x<-seq(from=0, to=1, by=0.01)
  K<-1
  chi<-rep$chi
  beta<-rep$beta
  
  uni_id<-unique(data$names)
  new_list<-list()
  for(i in 1:length(uni_id)){
    one_dex<-subset(data, names==uni_id[i])
    new_list[[i]]<-aggregate(one_dex$pa, by=list(one_dex$year), FUN=mean, na.rm=TRUE)
  }
  
  multi_inner <- Reduce(
    function(x, y, ...) plyr::join(x, y, by="Group.1", ...), 
    new_list
  )
  
  cols = c("black", "red", "purple")
  
  type3<-(K*new_x^beta[1])/((chi[1]^beta[1])+(new_x^beta[1]))
  plot(new_x, type3, type="l", lwd=2, ylab="Consumption Rate", xlab="Prey Density", las=1, ylim=c(0,1), cex.axis=1.25, cex.lab=1.5)
  if(length(rep$beta)>1){
    for(i in 2:length(rep$beta)){
      type3<-(K*new_x^beta[i])/((chi[i]^beta[i])+(new_x^beta[i]))
      lines(new_x, type3,lwd=2, col=cols[i])
    }
  }
  points(multi_inner[,3]~multi_inner[,2], pch=19)
  if(length(new_list)==3){points(multi_inner[,4]~multi_inner[,2], col="red")}
  if(length(new_list)==4){points(multi_inner[,5]~multi_inner[,2], col="purple")}
}



