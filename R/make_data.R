


#' Create basic dataframe for the model
#'
#'Combine presence/absence data from bottom trawl surveys, stomach contents, and an index for the different indices
#'
#' @param pa Vector of presence/absence data from all surveys (binary)
#' @param year Vector of years associated with each observation (numeric)
#' @param n Vector of the number of observations for each survey (integer)
#' @param id Vector of species identifiers for each data source (integer)
#' @param names Vector of names for each type of data, used for plotting (character)
#'
#' @return dataframe
#' @export
#'
#' @examples
#' make_data(pa=rbinom(n=20, size=1, prob=0.5), year=rep(seq(from=1, to=10, by=1),2), n=c(10,10), id =c(0,1,1), names=c("Trawl", "Call", "Full"))
make_data<-function(pa, year, n, id, names){
  fishy_dat<-data.frame(pa = pa, year = year, idex=rep(id, n), names = rep(names, n))
  return(fishy_dat)}