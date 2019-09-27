################################################################################
#' @title Fitting a pharmacokinetic model to concentration data
#' @description This function fits a pharmacokinetic model (dissolution profile) to time-concentration data using non-linear least squares regression.
#' The model can be chosen from a candidate set containing a First order, Hixson-Crowell,Higuchi, Weibull and a logistic model.
#' See Moellenhoff et al. (2018) <doi:10.1002/sim.7689> for details.
#'
#' @name fit_pharm_mod
#' @export
#' @param time a vector containing the time points of measurements
#' @param conc data frame or matrix containing the concentrations (see the example)
#' @param m model type. Built-in models are "firstorder",  "hixson",  "higuchi", "weibull" and "logistic"
#' @param plot plot of the model, default is TRUE.
#' @return A list containing the model type and the obtained parameters, further the RSS for all possible models. Furthermore a plot is given.
#' @examples
#' data(example_data)
#' conc1 <- select(filter(example_data,Group=="1"),-Tablet,-Group)
#' time <- c(10,15,20,30,45,60)
#' fit_pharm_mod(time,conc1,m="logistic")
#' @references Moellenhoff et al. (2018) <doi:10.1002/sim.7689>
################################################################################

fit_pharm_mod <- function(time,conc,m,plot=TRUE){
  #wrong specified data
  if (is.null(time)){
    stop("The time is not specified correctly")
  }
  if (is.null(conc)){
    stop("The concentration data is not specified correctly")
  }


  builtIn <- c("firstorder", "hixson", "higuchi", "weibull",
               "logistic");
  m1index <- match(m, builtIn);

  #wrong specification of the model
  if (is.na(m1index)) {
    stop("Invalid pharmacokinetic model specified")
  };

  p <- length(time)

  maxt <- max(time)
  mint <- min(time)

  n <- dim(conc)[1] #save sample size

  #preparations for fitting the models#
  x <- matrix(rep(time,n),ncol=p,byrow=TRUE)

  y <- as.matrix(conc,byrow=TRUE,dimnames=NULL)

  ##Objective functions##
  leastsquares_firstorder <- function(x,y){
    function(v){
      alpha1=v[1]
      sum((as.vector(y)-100*(1-exp(-alpha1*as.vector(x))))^2)}
  }
  leastsquares_hixson <- function(x,y){
    function(v){
      alpha1=v[1]
      sum((as.vector(y)-100*(1-(1-alpha1*as.vector(x)/4.6416)^3))^2)}
  }
  leastsquares_higuchi <- function(x,y){
    function(v){
      alpha1=v[1]
      sum((as.vector(y)-alpha1*as.vector(x)^0.5)^2)}
  }
  leastsquares_weibull <- function(x,y){
    function(v){
      k=v[1];beta=v[2]
      sum((as.vector(y)-100*(1-exp(-(as.vector(x)/k)^beta)))^2)}
  }
  leastsquares_logistic <- function(x,y){
    function(v){
      alpha=v[1];beta=v[2]
      sum((as.vector(y)-100*(exp(alpha+beta*log(as.vector(x))))/(1+exp(alpha+beta*log(as.vector(x)))))^2)}
  }

  #Matrix storing RSS
  RSS <- matrix(NA,nrow=5,ncol=1,dimnames=list(c("First order","Hixson-Crowell","Higuchi","Weibull","Logistic"),NULL))
  RSS[1,1] <- optimize(leastsquares_firstorder(x,y),c(0,1))$objective
  RSS[2,1] <- optimize(leastsquares_hixson(x,y),c(0,1))$objective
  RSS[3,1] <- optimize(leastsquares_higuchi(x,y),c(0,100))$objective
  RSS[4,1] <- optim(par=c(10,1),leastsquares_weibull(x,y))$value
  RSS[5,1] <- optim(par=c(-4,2),leastsquares_logistic(x,y))$value

  ###determining the model parameters###
  ## Model 1 ##
  if(m1index==1){
    leastsquares1 <- leastsquares_firstorder(x,y)
    coef1 <- optimize(leastsquares1,c(0,1))$minimum
    model.type <- "First order"
  }else if(m1index==2){
    leastsquares1 <- leastsquares_hixson(x,y)
    coef1 <- optimize(leastsquares1,c(0,1))$minimum
    model.type <- "Hixson-Crowell"
  }else if(m1index==3){
    leastsquares1 <- leastsquares_higuchi(x,y)
    coef1 <- optimize(leastsquares1,c(0,100))$minimum
    model.type <- "Higuchi"
  }else if(m1index==4){
    leastsquares1 <- leastsquares_weibull(x,y)
    coef1 <- optim(par=c(10,1),leastsquares1)$par
    model.type <- "Weibull"
  }else{
    leastsquares1 <- leastsquares_logistic(x,y)
    coef1 <- optim(par=c(-4,2),leastsquares1)$par
    model.type <- "Logistic"
  }

  ##determining the curves##
  #m1#
  if(m1index==1){
    modell1 <- function(d){100*(1-exp(-coef1[1]*d))}
  }else if(m1index==2){
    modell1 <- function(d){100*(1-(1-coef1[1]*d/4.6416)^3)}
  }else if(m1index==3){
    modell1 <- function(d){coef1[1]*d^0.5}
  }else if(m1index==4){
    modell1 <- function(d){100*(1-exp(-(d/coef1[1])^coef1[2]))}
  }else if(m1index==5){
    modell1 <- function(d){100*(exp(coef1[1]+coef1[2]*log(d))/(1+exp(coef1[1]+coef1[2]*log(d))))}
  }
  if(plot==TRUE) {
  curve(modell1,xlim=c(0,maxt),ylim=c(0,max(max(y),modell1(maxt))),xlab="Time",ylab="Percentage dissolved",main="Dissolution Curve")
  points(x,y)
  }
  #legend("bottomright",c(expression(m[1]),expression(m[2])),lty=c(1:3))
  return(list(type.of.model=model.type,estimated.model.param.=coef1, RSS=RSS))
}
