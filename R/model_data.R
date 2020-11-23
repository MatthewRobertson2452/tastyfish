#' Create data list for general functional response model
#'
#' @param fishy_dat Object output from the make_data() function
#' @param k The K parameter value for the general functional response model. Represents the upper asymptote of the curve.
#' @param chi The chi parameter value for the general functional response model. Represents the midpoint of the curve.
#' @param id An ID vector the length of the number of datasets used to represent different functional responses
#' @param n A vector the length of the number of datasets used where each number indicates the number of observations for that dataset
#'
#' @return A list of data to be used when generating the TMB model. The list includes the total number of observations (n), the number of years (nyrs),
#' the number of datasets (ndex), an index vector for the year of each observation (iyear), an index vector for the dataset associated with each observation (idex),
#' the k paramater, the chi parameter, a vector for all presence/absence data (pa), and an index vector identifying which functional response data will be associated with (t2dex)
#' @export
#'
#' @examples model_data(fishy_dat=fishy_dat, k=0.95, chi=0.45, id=c(0,1,1), n=c(length(df$pa),length(df1$pa),length(df2$pa)))
model_data<-function(fishy_dat, k, chi, id, n){
  tmb.data<-list(
    n=length(fishy_dat$pa),
    nyrs=length(unique(fishy_dat$year)),
    ndex=length(unique(fishy_dat$idex)),
    iyear=fishy_dat$year-min(fishy_dat$year),
    idex=fishy_dat$idex,
    k=k,
    chi=chi,
    pa=fishy_dat$pa, 
    t2dex=rep(id, times=n))
  if(min(tmb.data$t2dex)!=0){tmb.data$t2dex<-tmb.data$t2dex-min(tmb.data$t2dex)}
  return(tmb.data)}