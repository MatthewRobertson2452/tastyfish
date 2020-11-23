

#' Plot time series of index and all dataseries
#'
#' @param rep the report object from the TMB model
#' @param sdrep the variance report object from the TMB model
#' @param year a vector of the years used in the model
#' @param dat_names a vector of the names of the different data series used for the legend
#' @param ylim vector for y-limit in the plot (defaults to c(0,1))
#'
#' @return
#' @export
#'
#' @examples
plot_ts<-function(rep, sdrep, year, dat_names, ylim=c(0,1)){
  
  sec_mu_true<-rep$new_mu[,2]
  sec_mu_true[sec_mu_true <0.01] <- NA
  
  if(dim(rep$mu)[2]==3){
  third_mu_true<-rep$new_mu[,3]
  third_mu_true[third_mu_true <0.01] <- NA
  }
  
  low_rho = sdrep$value - sdrep$sd
  high_rho = sdrep$value + sdrep$sd
  
  plot(rep$iye~year, type="l", xlab="Year", ylab="Prey Index", ylim=ylim, lwd=2, cex.axis=1.25, cex.lab=1.5)
  polygon(c(year,rev(year)),c(low_rho,rev(high_rho)),border=NA,col="lightgrey")
  lines(rep$iye~year, lwd=2)
  lines(rep$new_mu[,1]~year, col="blue", lty=2)
  lines(sec_mu_true~year, col="red",lty=2)
  if(dim(rep$mu)[2]==3){lines(third_mu_true~year, col="purple",lty=2)}
  
  if(dim(rep$mu)[2]==2){col=c("black","blue","red")}
  if(dim(rep$mu)[2]==3){col=c("black","blue","red", "purple")}
  legend("topleft", legend=dat_names, 
         col=col, lwd=2, 
         lty=c(1,rep(2, length(dat_names))), bty="n")
}
