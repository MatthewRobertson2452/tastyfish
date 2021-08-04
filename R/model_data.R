#' Create data list for general functional response model
#'
#' @param fishy_dat Object output from the make_data() function
#' @param n A vector the length of the number of datasets used where each number indicates the number of observations for that dataset
#' @param type A character identifying whether you are creating data for a linear or nonlinear functional response model. Input "linear" for a linear model and "nonlinear" for a nonlinear model
#'
#' @return A list of data to be used when generating the TMB model. The list includes the total number of observations (n), the number of years (nyrs),
#' the number of datasets (ndex), an index vector for the year of each observation (iyear), an index vector for the dataset associated with each observation (idex), and
#' a vector for all presence/absence data (pa)
#' @export
#'
#' @examples model_data(fishy_dat=fishy_dat, n=c(length(df$pa),length(df1$pa),length(df2$pa)), type="nonlinear)
model_data<-function(fishy_dat, n, type){
  
  if(type=="linear"){
    tmb.data<-list(
      n=length(fishy_dat$pa),
      nyrs=length(unique(fishy_dat$year)),
      ndex=length(unique(fishy_dat$idex)),
      iyear=fishy_dat$year-min(fishy_dat$year),
      idex=fishy_dat$idex,
      pa=fishy_dat$pa)
  }
  
  if(type=="nonlinear"){
  tmb.data<-list(
    n=length(fishy_dat$pa),
    nyrs=length(unique(fishy_dat$year)),
    ndex=length(unique(fishy_dat$idex)),
    iyear=fishy_dat$year-min(fishy_dat$year),
    idex=fishy_dat$idex,
    k=1,
    pa=fishy_dat$pa)
  }
  
  return(tmb.data)}