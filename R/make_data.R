


#' Create basic dataframe for the model
#'
#'Combine presence/absence data from bottom trawl surveys, stomach contents, and an index for the different indices
#'
#' @param pa Vector of presence/absence data from all surveys (binary)
#' @param year Vector of years associated with each observation (numeric)
#' @param n Vector of the number of observations for each survey (integer)
#'
#' @return dataframe
#' @export
#'
#' @examples
#' make_data(pa=rbinom(n=20, size=1, prob=0.5), year=rep(seq(from=1, to=10, by=1),2), n=c(10,10))
make_data<-function(pa, year, n){
  fishy_dat<-data.frame(pa = pa, year = year, idex=rep(seq(from=1, to=length(n)), times=n)-1)
  return(fishy_dat)}