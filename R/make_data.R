


#' Create basic dataframe for the model
#'
#'Combine presence/absence data from bottom trawl surveys, stomach contents, and an index for the different indices
#'
#' @param pa Vector of presence/absence data from all surveys (binary)
#' @param year Vector of years associated with each observation (numeric)
#' @param nid Number of unique types of data used (integer)
#'
#' @return dataframe
#' @export
#'
#' @examples
#' make_data(pa=rbinom(n=20, size=1, prob=0.5), year=seq(from=1, to=10, by=1), nid=2)
make_data<-function(pa, year, nid){
  fishy_dat<-data.frame(pa = pa, year = rep(year, nid), idex=rep(seq(from=1, to=nid), each=length(year))-1)
  return(fishy_dat)}
