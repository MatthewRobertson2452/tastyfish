

#' Plot time series of index and all dataseries
#'
#' @param rep The report object from the TMB model
#' @param sdrep The variance report object from the TMB model
#' @param year A vector of the years used in the model
#' @param dat_names A vector of the names of the different data series used for the legend
#' @param ylim A vector for y-limit in the plot (defaults to c(0,1))
#'
#' @return plot
#' @export
#'
#' @examples plot_ts(rep=rep, sdrep=sdrep, year=seq(from=1995, to=2018, by=1), dat_names=c("Trawl","Full Stomach Contents","Called Stomach Contents"), ylim=c(0.2,0.6))
plot_ts<-function(tmb_data, rep, sdrep, year, dat_names, ylim=c(0,1)){
  
  dat_names<-c("Index", dat_names)
  
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
  
  multi_inner$firstdat<-((multi_inner[,3]*(rep$chi^rep$beta))/(1-multi_inner[,3]))^(1/rep$beta)
  if(length(new_list)==3){multi_inner$seconddat<-((multi_inner[,4]*(rep$chi^rep$beta))/(1-multi_inner[,4]))^(1/rep$beta)}
  

  low_rho = sdrep$value - sdrep$sd
  high_rho = sdrep$value + sdrep$sd

  multi_inner[,2]<-log(multi_inner[,2]/(1-multi_inner[,2]))
  multi_inner$firstdat<-log(multi_inner$firstdat/(1-multi_inner$firstdat))
  if(length(new_list)==3){multi_inner$seconddat<-log(multi_inner$seconddat/(1-multi_inner$seconddat))}
  
  n_dfs<-length(multi_inner)
  
  min_dat<-min(c(rep$ny, multi_inner[,2], multi_inner$firstdat, multi_inner$seconddat, low_rho, high_rho),na.rm=TRUE)-0.1
  max_dat<-max(c(rep$ny, multi_inner[,2], multi_inner$firstdat, multi_inner$seconddat, low_rho, high_rho),na.rm=TRUE)+0.1
    
  plot(rep$ny~year, type="l", xlab="Year", ylab="Prey Index", ylim=c(min_dat,max_dat), lwd=2, cex.axis=1.25, cex.lab=1.5)
  polygon(c(year,rev(year)),c(low_rho,rev(high_rho)),border=NA,col="lightgrey")
  lines(rep$ny~year, lwd=2)
  lines(multi_inner[,2]~year, col="blue", lty=2)
  lines(multi_inner$firstdat~year, col="red",lty=2)
  if(length(new_list)==3){lines(multi_inner$seconddat~year, col="purple",lty=2)}
  
  if(length(new_list)==2){col=c("black","blue","red")}
  if(length(new_list)==3){col=c("black","blue","red", "purple")}
  legend("topleft", legend=dat_names, 
         col=col, lwd=2, 
         lty=c(1,rep(2, length(dat_names))), bty="n")
}
