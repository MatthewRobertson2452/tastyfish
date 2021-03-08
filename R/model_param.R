#' Create parameter list for for general functional response model
#'
#' @param tmb_data TMB data object from the model_data() function
#' @param lbeta A logged starting value for the beta parameter
#' @param lchi A logged starting value for the chi parameter
#' @param lk A placeholder parameter for the linear model that allows it to run, since otherwise there wouldn't be any fixed effects parameters. Don't change from default
#' @param type A character identifying whether you are creating data for a linear or nonlinear functional response model. Input "linear" for a linear model and "nonlinear" for a nonlinear model
#'
#' @return A list of parameters for the functional response model, including a vector the length of the number of years (ny), 
#' and a starting value for the beta and chi parameter's
#' @export
#'
#' @examples model_param(tmb_data=tmb_data, lbeta=log(2), lchi=log(0.5), type="nonlinear")
model_param<-function(tmb_data, lbeta, lchi, type){
  
  if(type=="linear"){
    tmb_params<-list(
      ny = rep(0, tmb_data$nyrs),
      lk=0
    )}
  
  if(type=="nonlinear"){
  tmb_params<-list(
    ny = rep(0, tmb_data$nyrs),
    lbeta=lbeta,
    lchi=lchi
  )}
  return(tmb_params)}