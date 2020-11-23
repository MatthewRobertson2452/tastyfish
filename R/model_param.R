#' Create parameter list for TMB
#'
#' @param tmb_data 
#' @param lbeta 
#'
#' @return
#' @export
#'
#' @examples
model_param<-function(tmb_data, lbeta){
  tmb_params<-list(
    yr_tau_mat = matrix(0, nrow=tmb_data$nyrs, ncol=tmb_data$ndex),
    iye = rep(0, tmb_data$nyrs),
    lbeta=lbeta
  )
  return(tmb_params)}