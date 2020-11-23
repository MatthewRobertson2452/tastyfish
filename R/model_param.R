#' Create parameter list for for general functional response model
#'
#' @param tmb_data TMB data object from the model_data() function
#' @param lbeta A logged starting value for the beta paramater
#'
#' @return A list of parameters for the functional response model, including a matrix with dimensions by number of years and datatypes(yr_tau_mat), 
#' a vector the length of the number of years (iye), and a starting value for the beta parameter
#' @export
#'
#' @examples model_param(tmb_data=tmb_data, lbeta=log(5))
model_param<-function(tmb_data, lbeta){
  tmb_params<-list(
    yr_tau_mat = matrix(0, nrow=tmb_data$nyrs, ncol=tmb_data$ndex),
    iye = rep(0, tmb_data$nyrs),
    lbeta=lbeta
  )
  return(tmb_params)}