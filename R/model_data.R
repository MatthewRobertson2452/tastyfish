#' Create data list for TMB
#'
#' @param fishy_dat 
#' @param k 
#' @param chi 
#' @param id 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
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