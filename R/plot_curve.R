

#' Plot the functional response curve and associated data
#'
#' @param rep The report object from the TMB model
#' @param tmb_data The datalist used for the TMB model
#'
#' @return plot
#' @export
#'
#' @examples plot_curve(rep=rep, tmb_data=tmb_data)
plot_curve<-function(rep, tmb_data){
  new_x<-seq(from=0, to=1, by=0.01)
  K<-tmb_data$k
  chi<-rep$chi
  beta<-rep$beta
  type3<-(K*new_x^beta)/((chi^beta)+(new_x^beta))
  
  id<-unique(tmb_data$idex)
  tmb_df<-data.frame(id=tmb_data$idex, pa=tmb_data$pa, yr=tmb_data$iyear)
  new_list<-list()
  for(i in 1:length(id)){
    one_dex<-subset(tmb_df, id==i-1)
    new_list[[i]]<-aggregate(one_dex$pa, by=list(one_dex$yr), FUN=mean, na.rm=TRUE)
  }
  
  multi_inner <- Reduce(
    function(x, y, ...) plyr::join(x, y, by="Group.1", ...), 
    new_list
  )
  
  plot(new_x, type3, type="l", lwd=2, ylab="Consumption Rate", xlab="Prey Density", las=1, ylim=c(0,1), cex.axis=1.25, cex.lab=1.5)
  points(multi_inner[,3]~multi_inner[,2], pch=19)
  if(length(new_list)==3){points(multi_inner[,4]~multi_inner[,2], col="red")}
}



