#' Create parameter list for for general functional response model
#'
#' @param tmb_data TMB data object from the model_data() function
#' @param lbeta A logged starting value for the beta paramater
#'
#' @return A list of parameters for the functional response model, including a vector the length of the number of years (ny), 
#' and a starting value for the beta and chi parameter's
#' @export
#'
#' @examples model_param(tmb_data=tmb_data, lbeta=log(2), lchi=log(0.5))
model_param<-function(tmb_data, lbeta, lchi){
  tmb_params<-list(
    ny = rep(0, tmb_data$nyrs),
    lbeta=lbeta,
    lchi=lchi
  )
  return(tmb_params)}