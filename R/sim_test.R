################################################################################
#' @title Bootstrap test for the assessment of similarity of drug dissolutions profiles via maximum deviation
#'
#' @description Function for testing whether two dissolution profiles are similar concerning the
#' hypotheses \eqn{H_0: \max_{t\in\mathcal{T}} |m_1(t,\beta_1)-m_2(t,\beta_2)|\geq \epsilon\ vs.\
#' H_1: \max_{t\in\mathcal{T}} |m_1(t,\beta_1)-m_2(t,\beta_2)|< \epsilon.}
#'
#' $m_1$ and $m_2$ are pharmacokinetic models chosen from a candidate set containing a First order, Hixson-Crowell,Higuchi, Weibull and a logistic model.
#'
#' See Moellenhoff et al. (2018) <doi:10.1002/sim.7689> for details.
#'
#' @name sim_test
#' @export
#' @import alabama dplyr mvtnorm
#' @importFrom graphics curve legend points
#' @importFrom stats optim optimize ecdf
#' @param time1,time2 vectors containing the time points of measurements for each of the two formulations; if not further specified the time points are identical in both groups
#' @param conc1,conc2 data frames or matrices containing the concentrations obtained for each of the two formulations (see the example)
#' @param m1,m2 model types. Built-in models are "firstorder",  "hixson",  "higuchi", "weibull" and "logistic"
#' @param epsilon positive argument specifying the equivalence threshold (in \%), default is 10\% corresponding to an f2 of 50 according to current guidelines
#' @param B number of bootstrap replications. If missing, default value of B is 1000
#' @param plot if TRUE, a plot of the absolute difference curve of the two estimated models will be given. The default is FALSE.
#' @return A list containing the p.value, the types of models, the f2, the maximum absolute difference of the models, the estimated model parameters, the number of bootstrap replications and a summary of the bootstrap test statistic. Furthermore plots of the two models are given.
#' @examples
#' data(example_data)
#' conc1 <- select(filter(example_data,Group=="1"),-Tablet,-Group)
#' conc2 <- select(filter(example_data,Group=="2"),-Tablet,-Group)
#' time <- c(10,15,20,30,45,60)
#' sim_test(time1=time,time2=time,conc1=conc1,conc2=conc2,m1="logistic",m2="logistic",B=500,plot=TRUE)
#' @references Moellenhoff et al. (2018) <doi:10.1002/sim.7689>
#' @references EMA (2010) <https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-bioequivalence-rev1_en.pdf>
################################################################################

sim_test <- function(time1,time2=time1,conc1,conc2,m1,m2,epsilon=10,B=1000,plot=FALSE){
  #wrong specified data
  if (is.null(time1)|is.null(time2)){
    stop("The time is not specified correctly")
  }
  if (is.null(conc1)|is.null(conc2)){
    stop("The concentration data is not specified correctly")
  }

  p1 <- length(time1)
  p2 <- length(time2)

  maxt <- max(max(time1),max(time2))
  mint <- min(min(time1),min(time2))

  y_grid <- seq(mint,maxt,0.1) #creating a partition of the time range for searching the maximum

  n1 <- dim(conc1)[1] #save sample sizes
  n2 <- dim(conc2)[1]

  #preparations for fitting the models#
  x1 <- matrix(rep(time1,n1),ncol=p1,byrow=TRUE)
  x2 <- matrix(rep(time2,n2),ncol=p2,byrow=TRUE)

  y1 <- as.matrix(conc1,byrow=TRUE,dimnames=NULL)
  y2 <- as.matrix(conc2,byrow=TRUE,dimnames=NULL)

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


  builtIn <- c("firstorder", "hixson", "higuchi", "weibull",
               "logistic");
  m1index <- match(m1, builtIn);
  m2index <- match(m2, builtIn);

  #wrong specification of the model
  if (is.na(m1index) | is.na(m2index)) {
    stop("Invalid pharmacokinetic model specified")
  };

  ###determining the model parameters###
  ## Model 1 ##
  if(m1index==1){
    leastsquares1 <- leastsquares_firstorder(x1,y1)
    coef1 <- optimize(leastsquares1,c(0,1))$minimum
    model.type1 <- "First order"
  }else if(m1index==2){
    leastsquares1 <- leastsquares_hixson(x1,y1)
    coef1 <- optimize(leastsquares1,c(0,1))$minimum
    model.type1 <- "Hixson-Crowell"
  }else if(m1index==3){
    leastsquares1 <- leastsquares_higuchi(x1,y1)
    coef1 <- optimize(leastsquares1,c(0,100))$minimum
    model.type1 <- "Higuchi"
  }else if(m1index==4){
    leastsquares1 <- leastsquares_weibull(x1,y1)
    coef1 <- optim(par=c(10,1),leastsquares1)$par
    model.type1 <- "Weibull"
  }else{
    leastsquares1 <- leastsquares_logistic(x1,y1)
    coef1 <- optim(par=c(-4,2),leastsquares1)$par
    model.type1 <- "Logistic"
  }

  ##Model 2##
  if(m2index==1){
    leastsquares2 <- leastsquares_firstorder(x2,y2)
    coef2 <- optimize(leastsquares2,c(0,1))$minimum
    model.type2 <- "First order"
  }else if(m2index==2){
    leastsquares2 <- leastsquares_hixson(x2,y2)
    coef2 <- optimize(leastsquares2,c(0,1))$minimum
    model.type2 <- "Hixson-Crowell"
  }else if(m2index==3){
    leastsquares2 <- leastsquares_higuchi(x2,y2)
    coef2 <- optimize(leastsquares2,c(0,100))$minimum
    model.type2 <- "Higuchi"
  }else if(m2index==4){
    leastsquares2=leastsquares_weibull(x2,y2)
    coef2 <- optim(par=c(10,1),leastsquares2)$par
    model.type2 <- "Weibull"
  }else{
    leastsquares2=leastsquares_logistic(x2,y2)
    coef2 <- optim(par=c(-4,2),leastsquares2)$par
    model.type2 <- "Logistic"
  }

  ##determining the curves##
  #m1#
  if(m1index==1){
    model1 <- function(v){function(d){100*(1-exp(-v[1]*d))}}
  }else if(m1index==2){
    model1 <- function(v){function(d){100*(1-(1-v[1]*d/4.6416)^3)}}
  }else if(m1index==3){
    model1 <- function(v){function(d){v[1]*d^0.5}}
  }else if(m1index==4){
    model1 <- function(v){function(d){100*(1-exp(-(d/v[1])^v[2]))}}
  }else if(m1index==5){
    model1 <- function(v){function(d){100*(exp(v[1]+v[2]*log(d))/(1+exp(v[1]+v[2]*log(d))))}}
  }
  curve1 <- model1(coef1)

  ##m2##
  if(m2index==1){
    model2 <- function(v){function(d){100*(1-exp(-v[1]*d))}}
  }else if(m2index==2){
    model2 <- function(v){function(d){100*(1-(1-v[1]*d/4.6416)^3)}}
  }else if(m2index==3){
    model2 <- function(v){function(d){v[1]*d^0.5}}
  }else if(m2index==4){
    model2 <- function(v){function(d){100*(1-exp(-(d/v[1])^v[2]))}}
  }else if(m2index==5){
    model2 <- function(v){function(d){100*(exp(v[1]+v[2]*log(d))/(1+exp(v[1]+v[2]*log(d))))}}
  }
  curve2 <- model2(coef2)

  #plot of the data and the curves
  curve(curve1,xlim=c(0,maxt),ylim=c(0,max(max(max(curve1(y_grid)),max(curve2(y_grid)),max(y1),max(y2)))),xlab="Time",ylab="Percentage dissolved",main="Dissolution Curves")
  curve(curve2,xlim=c(0,maxt),add=TRUE,lty=3)
  points(x1,y1)
  points(x2,y2,pch=3)
  legend("bottomright",c(expression(m[1]),expression(m[2])),lty=c(1:3))

  #small B
  if (B<200) {
    warning("Warning: A larger B should be choosen for higher accuracy of the test.")
  }

  #calculate the test statistic
  difference <- function(d){abs(curve1(d)-curve2(d))}
  t.stat <- max(difference(y_grid))

  #obtain estimates of the covariance matrices of the data
  eps_1 <- matrix(as.vector(t(y1))-curve1(rep(time1,n1)),nrow=n1,ncol=p1,byrow=TRUE) #Matrix mit den Fehlervekoren (Anzahl Zeilen=Anzahl Patienten)
  epsquer <- colMeans(eps_1)
  v <- matrix(NA,nrow=n1,ncol=p1)
  vt <- list()
  for (i in 1:n1){
    v[i,] <- eps_1[i,]-epsquer
  }
  for (i in 1:n1){
    vt[i] <- list(v[i,]%*%t(v[i,]))
  }
  S1hat <- 1/(n1-2)*Reduce('+', vt)

  eps_2 <- matrix(as.vector(t(y2))-curve2(rep(time2,n2)),nrow=n2,ncol=p2,byrow=TRUE) #Matrix mit den Fehlervekoren (Anzahl Zeilen=Anzahl Patienten)
  epsquer2 <- colMeans(eps_2)
  v2 <- matrix(NA,nrow=n2,ncol=p2)
  vt2 <- list()
  for (i in 1:n2){
    v2[i,] <- eps_2[i,]-epsquer2
  }
  for (i in 1:n2){
    vt2[i] <- list(v2[i,]%*%t(v2[i,]))
  }
  S2hat <- 1/(n2-2)*Reduce('+', vt2)

  nop1 <- length(coef1) #length of parameters
  nop2 <- length(coef2)

  ### Bootstrap Test ###
  #constrained optimization
  if (t.stat>=epsilon){
    minimum <- c(coef1,coef2)} else {
      leastsquares <- function(v){
        alpha <- v[1:nop1]; beta <- v[(nop1+1):length(v)];
        leastsquares1(alpha)+leastsquares2(beta)
      }
      softmax <- function(v){
        alpha <- v[1:nop1]; beta <- v[(nop1+1):length(v)];
        dff <- function(x){abs(model1(alpha)(x)-model2(beta)(x))};
        sum(dff(y_grid)*exp(50*dff(y_grid)))/(sum(exp(50*dff(y_grid))))-epsilon
      }
      minimum <- auglag(par=c(coef1,coef2),leastsquares,heq=softmax,control.outer=list(trace=0))$par
    }

  #bootstrap
  boot <- vector()
  for (i in 1:B)
  {
    y1star <- model1(minimum[1:nop1])(x1)+rmvnorm(n1,rep(0,p1),sigma=S1hat)
    y2star <- model2(minimum[(nop1+1):length(minimum)])(x2)+rmvnorm(n2,rep(0,p2),sigma=S2hat)
    #fitting the bootstrap models
    ###determining the model parameters###
    ## Model 1 ##
    coef1star <- fit_pharm_mod(time=time1,conc=y1star,m=m1,plot=FALSE)$estimated.model.param
    ##Model 2##
    coef2star <- fit_pharm_mod(time=time2,conc=y2star,m=m2,plot=FALSE)$estimated.model.param
    boot[i] <- max(abs(model1(coef1star)(y_grid)-model2(coef2star)(y_grid)))#calculating the bootstrap test statistic
  }
  #calculating the p-value by using the ecdf of the bootstrap distribution
  fn <- ecdf(boot)
  pval <- fn(t.stat)
  if(plot==TRUE) {
    plot(y_grid,abs(model1(coef1)(y_grid)-model2(coef2)(y_grid)),type="l",xlim=c(0,maxt),xlab="Time",ylab="Percentage dissolved",main="Absolute difference curve of the two fitted models");
  }

  #f2
  if (n1<12|n2<12) {
    warning("Warning: f2 requires at least 12 samples.")
  }
  if (any(time1!=time2)) {
    warning("Warning: f2 can only be calculated if the sampling times are identical.")
  }
  if(all(time1==time2)){
  return(list(p.value=pval,types.of.models=c(model.type1,model.type2),max.abs.difference=t.stat,f2=f2(conc1,conc2),estimated.model.param.=c(coef1,coef2),bootstrap.replications=B,summary.boot=summary(boot)))
  } else {
    return(list(p.value=pval,types.of.models=c(model.type1,model.type2),max.abs.difference=t.stat,estimated.model.param.=c(coef1,coef2),bootstrap.replications=B,summary.boot=summary(boot)))
      }
}
