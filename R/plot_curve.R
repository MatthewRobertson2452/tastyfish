

#' Plot the functional response curve and associated data
#'
#' @param rep the report object from the TMB model
#' @param tmb_data the datalist used for the TMB model
#'
#' @return
#' @export
#'
#' @examples
plot_curve<-function(rep, tmb_data){
  new_x<-seq(from=0, to=1, by=0.01)
  K<-tmb_data$k
  chi<-tmb_data$chi
  beta<-rep$beta
  type3<-(K*new_x^beta)/((chi^beta)+(new_x^beta))
  plot(new_x, type3, type="l", lwd=2, ylab="Consumption Rate", xlab="Prey Density", las=1, ylim=c(0,1), cex.axis=1.25, cex.lab=1.5)
  points(rep$mu[,2]~rep$mu[,1], pch=19)
  if(dim(rep$mu)[2]==3){points(rep$mu[,3]~rep$mu[,1], col="red")}
}
