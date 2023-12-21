#############################################################################
#### Asymetry in lomax power distribution ###############################
###########################################################################

F_dlomax <- function(x){
  prob = vector(length = length(x))
  
  prob[which(x>0)] = 1 - 0.5*(1/(1+x[x>0]))
  prob[which(x<=0)] = 0.5*(1/(1 - x[x<=0]))
  
  return(prob)
}


Fplomax <- function(x, loglambda){
  lambda = exp(loglambda)
  p = F_dlomax(x)^lambda
  return(p)
}

Finvplomax <- function(x, loglambda){
  lambda = exp(loglambda)
  p = 1 - F_dlomax(-x)^lambda
  return(p)
}

Qlomax <- function(q,lambda){
  if(q>0.5){
    x = 1/(2*(1-q)) - 1 
  }else{
    x = 1 - 1/(2*q)
  }
  return(x)
}

Qplomax <- function(q,lambda){
  x = Qlomax(q^(1/lambda))
  return(x)
}

Qrplomax <- function(q,lambda){
  x = - Qlomax((1-q)^(1/lambda))
  return(x)
}

asymetry <- function(dist, lambda = 1){
  a = ((dist(0.875,lambda) - dist(0.5,lambda)) - (dist(0.5,lambda) - dist(0.125,lambda)))/(dist(0.875,lambda) - dist(0.125,lambda))
  return(a)
}


lambda <- seq(0,1, by = 0.001)

df <- data.frame('lambda' = lambda, 'dist' = 'PDLomax','assimetria' = unlist(lapply(lambda, function(x) asymetry(Qplomax,x))))
assimetria = rbind(df, data.frame('lambda' = lambda, 'dist' = 'RPDLomax','assimetria' = unlist(lapply(lambda, function(x) asymetry(Qrplomax,x))))) 

library(ggplot2)
ggplot(assimetria, aes(x=lambda, y=assimetria, group = dist)) +geom_line(aes(linetype = dist))+
  xlab("\u03BB")+ ylab(bquote(A[0]))+labs(linetype = "Distribuição")+
  theme_bw()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "Assimetria2.png", width = 6, height = 4)

lambda <- seq(0,10, by = 0.01)

df <- data.frame('lambda' = lambda, 'dist' = 'PDLomax','assimetria' = unlist(lapply(lambda, function(x) asymetry(Qplomax,x))))
assimetria = rbind(df, data.frame('lambda' = lambda, 'dist' = 'RPDLomax','assimetria' = unlist(lapply(lambda, function(x) asymetry(Qrplomax,x))))) 

library(ggplot2)
ggplot(assimetria, aes(x=lambda, y=assimetria, group = dist)) +geom_line(aes(linetype = dist))+
  xlab("\u03BB")+ ylab(bquote(A[0]))+labs(linetype = "Distribuição")+
  theme_bw()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "Assimetria1.png", width = 6, height = 4)

