##############################################################################
#### métricas de ajuste ######################################################
#############################################################################

## DIC
dic <- function(x, y, logLikelihood0, postSamples){
  
  logLikelihood <- function(theta) logLikelihood0(theta, y, x)
  
  thetaBayes <- colMeans(postSamples)
  
  dbar <- mean(-2 * apply(postSamples, 1, logLikelihood))
  dhat <- -2 * logLikelihood(thetaBayes)
  dic <- 2 * dbar - dhat
  
  measures = list("dbar"=dbar, "dhat"=dhat, "rho" = dbar-dhat,"dic"=dic)
  return(measures)
}

## EAIC
eaic <- function(x,y, logLikelihood0, postSamples){
  
  logLikelihood <- function(theta) logLikelihood0(theta, y, x)
  
  dbar <- mean(-2 * apply(postSamples, 1, logLikelihood) )
  eaic <- dbar + 2*ncol(postSamples)
  
  return(eaic)
}

## EBIC
ebic <- function(x,y, logLikelihood0, postSamples){
  n <- length(y)
  
  logLikelihood <- function(theta) logLikelihood0(theta, y, x)
  
  dbar <- mean(-2 * apply(postSamples, 1, logLikelihood) )
  ebic <- dbar + ncol(postSamples)*log(n)
  
  return(ebic)
}


##################################################################################
########## Distributions #########################################################
##################################################################################
#Double lomax
F_dlomax <- function(x){
  prob = vector(length = length(x))
  
  prob[which(x>0)] = 1 - 0.5*(1/(1+x[x>0]))
  prob[which(x<=0)] = 0.5*(1/(1 - x[x<=0]))
  
  return(prob)
}


loglike_dlomax <- function(par,y,x){
  eta <-x %*% par
  p <- F_dlomax(eta)
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero) 
}

#Power Lomax
Fplomax <- function(x, loglambda){
  lambda = exp(loglambda)
  p = F_dlomax(x)^lambda
  return(p)
}

loglike_dplomax <- function(par,y,x){
  eta <-x %*% par[-length(par)]
  p <- Fplomax(eta, par[length(par)])
  logvero = sum(y*log(p) + (1-y)*log(1-p))
  return(logvero) 
}

## Reverse Power Lomax
Finvplomax <- function(x, loglambda){
  lambda = exp(loglambda)
  p = 1 - F_dlomax(-x)^lambda
  return(p)
}

loglike_dpinvlomax <- function(par,y,x){
  eta <-x %*% par[-length(par)]
  p <- Finvplomax(eta, par[length(par)])
  logvero = sum(y*log(p) + (1-y)*log(1-p))
  return(logvero) 
}

#Logística
loglike_log <- function(par,y,x) {
  xbeta=x%*%par	
  p=exp(xbeta)/(1+exp(xbeta)) 	
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero)
}

#Probit
loglike_probit <- function(par,y,x) {
  xbeta=x%*%par	
  p=pnorm(xbeta)	
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero)
}

### Cloglog
loglike_cloglog <- function(par,y,x) {
  xbeta=x%*%par	
  p=1-exp(-exp(xbeta)) 	
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero)
}

### loglog
loglike_loglog <- function(par,y,x) {
  xbeta=x%*%par	
  p=exp(-exp(-xbeta)) 	
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero)
}


### Cauchit
loglike_cauchit <- function(par,y,x) {
  xbeta=x%*%par	
  p=pcauchy(xbeta)	
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero)
}

