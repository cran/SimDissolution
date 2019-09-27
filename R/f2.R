################################################################################
#' @title Function for computing the f2
#' @description Function for computing the f2, time points have to be identical.
#' Validity criteria of the f2 have to be be checked in advance.
#' See Moellenhoff et al. (2018) <doi:10.1002/sim.7689>
#'
#' @name f2
#' @export
#' @param conc1,conc2 data frames containing the concentrations obtained for each of the two formulations
#' @return a single value for the f2
#' @examples
#' data(example_data)
#' conc1<-select(filter(example_data,Group=="1"),-Tablet,-Group)
#' conc2<-select(filter(example_data,Group=="2"),-Tablet,-Group)
#' f2(conc1=conc1,conc2=conc2)
#' @references Moellenhoff et al. (2018) <doi:10.1002/sim.7689>
################################################################################

f2 <- function(conc1,conc2){
  mu1=colMeans(conc1)
  mu2=colMeans(conc2)
  d=mean((mu1-mu2)^2)
  f2=50*log(100/sqrt((1+d)))/log(10)
  return(f2)
}
