##############################################################################
##### Curvas Logísticas ######################################################
#################################################################################
require(ggplot2)
require(ggsave)

### link functions
logit <- function(eta){exp(eta)/(1+exp(eta))}
probit <- function(eta){pnorm(eta)}
cauchit <-  function(eta){0.5 + atan(eta)/pi}
cloglog <- function(eta){exp(-exp(-eta))}
loglog <- function(eta){1-exp(-exp(eta))}

### gráfico com as link functions
eta <- seq(-10,10, by = 0.1)

probs <- data.frame(eta = eta, Link = rep(c("Logit","Probit", "Cauchit", "Cloglog", "Loglog"), each = length(eta)),
                    probs = c(logit(eta), probit(eta), cauchit(eta), cloglog(eta), loglog(eta)))
  
theme_set(theme_bw())
ggplot(probs, aes(x=eta, y=probs, group = Link)) +geom_line(aes(linetype=Link))+
  xlab(expression(eta[i]))+ ylab(expression(p[i]))
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F1_ligacoes_comuns.png", width = 6, height = 4)

################################################################################
#### Lomax distribution ###################################################
#################################################################################

f_lomax <- function(x, lambda, alpha){
  prob <- alpha/lambda * (1+ x/lambda)^(-alpha-1)
  return(prob)
}

F_lomax <- function(x, lambda, alpha){
  cum_prob <- 1 - ( 1 + x/lambda)^(-alpha)
  return(cum_prob)
}

Fp_lomax <- function(x, lambda, alpha,gamma){
  cum_prob <- (1 - ( 1 + x/lambda)^(-alpha))^gamma
  return(cum_prob)
}

Finv_lomax <- function(p, lambda, alpha){
  eta <- lambda * ((1-p)^alpha - 1)
  return(eta)
}


x <- seq(0,10, by=0.01)
plot(F_lomax(x, 5, 4), type = 'l',x = x)


library(Renext)
plomax(1,2,2)
F_lomax(1,2,2)
dlomax(1,2,2)
####################################################################
##### Gráficos da Lomax ############################################
######################################################################
# Distribuição Lomax
# diferentes alpha
x <- seq(0.01,10, by = 0.1)
alpha1 = f_lomax(x,lambda=2,alpha=1)
alpha2 = f_lomax(x,lambda=2,alpha=2)
alpha3 = f_lomax(x,lambda=2,alpha=4)
alpha4 = f_lomax(x,lambda=2,alpha=6)

probs <- data.frame(x=x, Link = rep(c("\u03b1=1","\u03b1=2", "\u03b1=4","\u03b1=6"), each = length(x)),
                    probs = c(alpha1, alpha2, alpha3, alpha4))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDP")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F2_densidadeLomaxAlpha.png", width = 11, height = 7)

#Diferentes lambda
x <- seq(0.01,10, by = 0.1)
alpha1 = f_lomax(x,lambda=2,alpha=2)
alpha2 = f_lomax(x,lambda=4,alpha=2)
alpha3 = f_lomax(x,lambda=6,alpha=2)
alpha4 = f_lomax(x,lambda=8,alpha=2)

probs <- data.frame(x=x, Link = rep(c("\u03BB=2","\u03BB=4", "\u03BB=6","\u03BB=8"), each = length(x)),
                    probs = c(alpha1, alpha2, alpha3, alpha4))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDP")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F3_densidadeLomaxLambda.png", width = 11, height = 7)

#Distribuição Acumulada Lomax
# diferentes alpha
x <- seq(0.01,10, by = 0.1)
alpha1 = F_lomax(x,lambda=2,alpha=1)
alpha2 = F_lomax(x,lambda=2,alpha=2)
alpha3 = F_lomax(x,lambda=2,alpha=4)
alpha4 = F_lomax(x,lambda=2,alpha=6)

probs <- data.frame(x=x, Link = rep(c("\u03b1=1","\u03b1=2", "\u03b1=4","\u03b1=6"), each = length(x)),
                    probs = c(alpha1, alpha2, alpha3, alpha4))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDP")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F4_acumuladaLomaxAlpha.png", width = 11, height = 7)

#Diferentes lambda
x <- seq(0.01,10, by = 0.1)
alpha1 = F_lomax(x,lambda=2,alpha=2)
alpha2 = F_lomax(x,lambda=4,alpha=2)
alpha3 = F_lomax(x,lambda=6,alpha=2)
alpha4 = F_lomax(x,lambda=8,alpha=2)

probs <- data.frame(x=x, Link = rep(c("\u03BB=2","\u03BB=4", "\u03BB=6","\u03BB=8"), each = length(x)),
                    probs = c(alpha1, alpha2, alpha3, alpha4))
library(ggplot2)
ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDP")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F5_acumuladaLomaxLambda.png", width = 11, height = 7)


#########################################################################################
############ Power Lomax ##############################################################
##########################################################################################
fp_lomax <- function(x, lambda, alpha, gamma){
  prob <- gamma*(alpha/lambda * (1+ x/lambda)^(-alpha-1))*(1 - ( 1 + x/lambda)^(-alpha))^(gamma-1)
  return(prob)
}

Fp_lomax <- function(x, lambda, alpha, gamma){
  cum_prob <- (1 - ( 1 + x/lambda)^(-alpha))^gamma
  return(cum_prob)
}
########################################
######### PDF ########################333
x <- seq(0.01,20, by = 0.1)
gamma1 = fp_lomax(x,lambda=1,alpha=2, gamma =70)
gamma2 = fp_lomax(x,lambda=4,alpha=2, gamma =2)
gamma3 = fp_lomax(x,lambda=4,alpha=2, gamma =3)
gamma4 = fp_lomax(x,lambda=4,alpha=2, gamma =4)
gamma5 = fp_lomax(x,lambda=4,alpha=2, gamma =5)

probs <- data.frame(x=x, Link = rep(c("\u03B3=1","\u03B3=2", "\u03B3=3","\u03B3=4", "\u03B3=5"), each = length(x)),
                    probs = c(gamma1, gamma2, gamma3, gamma4, gamma5))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDP")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F6_densidadePLomaxLambda.png", width = 11, height = 7)

########################################
######### CDF ########################333
x <- seq(0.01,50, by = 0.1)
gamma1 = Fp_lomax(x,lambda=1,alpha=2, gamma =70)
gamma2 = Fp_lomax(x,lambda=4,alpha=2, gamma =2)
gamma3 = Fp_lomax(x,lambda=4,alpha=2, gamma =3)
gamma4 = Fp_lomax(x,lambda=4,alpha=2, gamma =4)
gamma5 = Fp_lomax(x,lambda=4,alpha=2, gamma =5)

probs <- data.frame(x=x, Link = rep(c("\u03B3=1","\u03B3=2", "\u03B3=3","\u03B3=4", "\u03B3=5"), each = length(x)),
                    probs = c(gamma1, gamma2, gamma3, gamma4, gamma5))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDA")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F7_acumuladaPLomaxLambda.png", width = 11, height = 7)


##############################################################
#### Máxima Verossimilhança ##################################
################################################################

#param = (alpha, lambda)
loglike_lomax <- function(param, x, y){
  eta <- x %*% as.matrix(param[1:2])
  mu <- 1 - (1 + (eta/param[3]))^(-param[4])
  loglike <- y*log(mu) + (1-y)*log(1-mu)
  return(sum(loglike))
}

########################################################################
#### Estudo de Simulação ###########################################3
#######################################################################
beta0 = 0 
beta1 = 1

n = 1000
r = 100
x  = runif(1000, 0, 20)
eta = beta0 + beta1 * x
alpha = 2
p <- 1 - (1 + eta)^(-alpha)
y = rbinom(1000,1,p)
dados <- cbind("x" = x,"y" = y)

loglike_lomax <- function(param, y, x){
  eta <-x%*%as.matrix(param)
  p <- 1 - (1 + eta)^(-2)
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero)
}


m.log = try(glm(formula = y ~ x, family = binomial(link = "logit"), data=data.frame(dados)), TRUE)

est.DL   = try(maxLik(logLik=loglike_lomax, start=c(m.log$coeff[1], m.log$coeff[2]), y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS"), TRUE) 

#############################
## Estimando alpha ##########
#############################
beta0 = 0
beta1 = 1

n = 1000
r = 100
x  = runif(1000, 0, 200)
eta = beta0 + beta1 * x
p <- 1 - (1 + eta)^(-0.75)
y = rbinom(1000,1,p)
dados <- cbind("x" = x,"y" = y)

loglike_lomax <- function(param, y, x){
  eta <-x%*%as.matrix(param[-length(param)])
  p <- 1 - (1 + eta)^(-param[length(param)])
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero) 
}


m.log = try(glm(formula = y ~ x, family = binomial(link = "logit"), data=data.frame(dados)), TRUE)

mom = 1+1/mean(y)

start = c(m.log$coefficients[1],m.log$coefficients[2],mom)

estimacao = maxLik(logLik=loglike_lomax,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS")
estimacao$estimate

############
# Simulação 
###########

r= 100
n=1000
parametros = data.frame(matrix(ncol=3,nrow=0))
colnames(parametros) = c("beta0", "beta1", "alpha")

for (i in 1:r){
    set.seed(i)
    x  = runif(n, 0,200)
    eta = beta0 + beta1 * x
    alpha = 0.5
    p <- 1 - ( 1 + eta)^(-alpha)
    y = rbinom(1000,1,p)
    dados <- data.frame("x" = x,"y" = y)
    
    m.log = try(glm(formula = y ~ x, family = binomial(link = "probit"), data=dados), TRUE)
    mom = 1+1/mean(y)
    start = c(m.log$coefficients[1],m.log$coefficients[2],mom)
    estim = maxLik(logLik=loglike_lomax,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS")
    parametros[i,] = estim$estimate
}

colMeans(parametros)

##################################333
### Estimação potência de Lomax #####
#####################################
beta0 = 2 
beta1 = 2

n = 1000
r = 100
x  = runif(1000, 0,20)
eta = beta0 + beta1 * x
alpha = 0.75
gamma = 2
p <- (1 - ( 1 + eta)^(-alpha))^gamma
y = rbinom(1000,1,p)
dados <- cbind("x" = x,"y" = y)



m.log = try(glm(formula = y ~ x, family = binomial(link = "logit"), data=data.frame(dados)), TRUE)

mom = 1+1/mean(y)

start = c(m.log$coefficients[1], m.log$coefficients[2],1,1)

maxLik(logLik=loglike_plomax,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="SANN")

############
# Simulação 
###########

r= 100
n=1000
parametros = data.frame(matrix(ncol=4,nrow=0))
colnames(parametros) = c("lambda","alpha","beta0", "beta1")

for (i in 1:r){
  set.seed(i)
  x  = runif(n, 0,200)
  eta = beta0 + beta1 * x
  alpha = 0.5
  gamma = 2
  p <- (1 - ( 1 + eta)^(-alpha))^gamma
  y = rbinom(1000,1,p)
  dados <- data.frame("x" = x,"y" = y)
  
  m.log = try(glm(formula = y ~ x, family = binomial(link = "cloglog"), data=dados), TRUE)
  #mom = 1+1/mean(y)
  start = c(m.log$coefficients[1],m.log$coefficients[2],1,1)
  estim = maxLik(logLik=loglike_plomax,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS")
  parametros[i,] = estim$estimate
}

colMeans(parametros)


##############################################################################
## Proporção 0 e 1 ####################################################3
############################################################################

n = 1000
r = 100
x  = runif(1000, 0,20)
eta = beta0 + beta1 * x
alpha = 0.75
gamma = 1
p <- (1 - ( 1 + eta)^(-alpha))^gamma
y = rbinom(1000,1,p)
table(y)

prop.table(table(y))

prop = data.frame(matrix(ncol=2, nrow=0))
colnames(prop) = c("prop", "gamma")
for (j in 1:4){  
  for (i in 1:r){
    set.seed(i)
    x  = runif(1000, 0,20)
    eta = beta0 + beta1 * x
    alpha = 0.75
    gamma = j
    p <- (1 - ( 1 + eta)^(-alpha))^gamma
    y = rbinom(1000,1,p)
    prop[i+r*(j-1),] = c(mean(y), j)
  }
}
prop$gamma = factor(prop$gamma)

require(ggplot2)
ggplot(prop, aes(y=prop, x=gamma)) + 
  geom_boxplot()+ylab('Proporções')+xlab('Potência')+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F2_densidadeLomaxAlpha.png", width = 11, height = 7)

#########################################################################################
### novos testes de distribuição ####################################################
######################################################################################
#Birnbaun Saunders
#nakagami 
#Hyperbolic secant distribution

## Birnbaun Sanders

F_bs <- function(x, alpha = 0.1, beta){
  cum_prob <- pnorm(1/alpha * (sqrt(x/beta) - sqrt(beta/x)))
  return(cum_prob)
}
########################################
######### PDF ########################333
gamma1 = F_bs(x, beta =0.5)
gamma2 = F_bs(x, alpha =0.3,beta =0.71)
gamma3 = F_bs(x, beta =1)
gamma4 = F_bs(x, beta =2)
gamma5 = F_bs(x, beta =3)

probs <- data.frame(x=x, Link = rep(c("\u03B3=1","\u03B3=2", "\u03B3=3","\u03B3=4", "\u03B3=5"), each = length(x)),
                    probs = c(gamma1, gamma2, gamma3, gamma4, gamma5))
library(ggplot2)
ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDP")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F6_densidadePLomaxLambda.png", width = 11, height = 7)

## Verossimilhança
loglike_bs <- function(param, y, x){
  eta <-as.matrix(x)%*%as.matrix(param[-length(param)])
  p <- F_bs(eta, beta= param[length(param)])
  logvero=sum(y*log(p) + (1-y)*log(1-p))
  return(logvero) 
}

### Simulação

beta0 = 30
beta1 = 20

n = 1000
x  = runif(1000, 0.5,100)
eta = beta0 + beta1 * x
beta = 450
p <- F_bs(eta, beta = beta)
y = rbinom(1000,1,p)
dados <- cbind("x" = x,"y" = y)


loglike_bs <- function(param, y, x){
  eta <-x%*%param[-length(param)]
  p <- F_bs(eta, beta= param[length(param)])
  logvero=log(prod(p^y*(1-p)^(1-y))+0.000000000001)
  return(logvero) 
}


logvero.L=function(param,y,x) {
  xbeta=x%*%as.matrix(param)	
  p=exp(xbeta)/(1+exp(xbeta)) 	
  logvero=log(prod(p^y*(1-p)^(1-y)))
  return(logvero)
}

start = c(1,1)
A <- diag(1,2)
B <- c(-0.1,-0.1)
m.log = maxLik(logLik=logvero.L,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS",
       constraints=list(ineqA=A, ineqB=B))



start = c(m.log$estimate[1],m.log$estimate[2],50)
A <- diag(1,3)
B <- c(-1e-9,-1e-9,-20)
maxLik(logLik=loglike_bs,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS",
       constraints=list(ineqA=A, ineqB=B))

##################################33
###  Double Weibul ################
####################################
eta <- seq(-10,10, by = 0.1)

c=2
F_dweibul <- function(x,c){
  prob = vector(length = length(x))
  
  prob[which(x>0)] = 1 - 0.5*exp(-abs(x[x>0])^c)
  prob[which(x<=0)] = 0.5*exp(-abs(x[x<=0])^c)
  
  return(prob)
}

plot(F_dweibul(eta,c =10))

r= 100
n=1000
parametros = data.frame(matrix(ncol=3,nrow=0))
colnames(parametros) = c("beta0", "beta1", "c")
beta0 = 0
beta1 = 1
c = 2

loglike_dweibul <- function(x, y, param){
  eta <-x %*% param[-length(param)]
  p <- F_dweibul(eta, c = param[length(param)])
  logvero = log(prod(p^y * (1-p)^(1-y)))
  return(logvero) 
}


for (i in 1:r){
  set.seed(i)
  x  = runif(n, -10, 10)
  eta = beta0 + beta1 * x
  p <- F_dweibul(eta,c =c)
  y = rbinom(1000,1,p)
  dados <- data.frame("x" = x,"y" = y)
  
  m.log = try(glm(formula = y ~ x, family = binomial(link = "logit"), data=dados), TRUE)

  start = c(m.log$coefficients[1],m.log$coefficients[2],1)
  estim = maxLik(logLik=loglike_dweibul,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS")
  parametros[i,] = estim$estimate
}

colMeans(parametros)

#################3#########################33
## Função de Densidade da Distribuição ###
##########################################
f_dweibul <- function(x, c){0.5*abs(x)^(c-1)*exp(-abs(x)^c)}

x <- seq(-10,10, by = 0.1)
c1 = f_dweibul(x, c=0.75)
c2 = f_dweibul(x, c=1)
c3 = f_dweibul(x,c =1.5)
c4 = f_dweibul(x,c =2)


probs <- data.frame(x=x, Link = rep(c("c=0.75","c=1","c=1.5", "c=2"), each = length(x)),
                    probs = c(c1, c2, c3,c4))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDA")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F10_densidadeDweibul.png", width = 11, height = 7)

################
#### Link Double Weibul
##############################

### link functions
probit <- function(eta){pnorm(eta)}

### gráfico com as link functions
eta <- seq(-15,15, by = 0.1)

probs <- data.frame(eta = eta, Link = rep(c("logit","Dweibul(c =4)", "Dweibul(c =0.5)"), each = length(eta)),
                    probs = c(logit(eta), F_dweibul(eta, c=4),  F_dweibul(eta, c=0.5)))

ggplot(probs, aes(x=eta, y=probs, group = Link)) +geom_line(aes(linetype=Link))+
  xlab(expression(eta[i]))+ ylab(expression(p[i]))+
  theme_bw()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F11_weibul_link_comparation.png", width = 11, height = 7)


####################################################################
### Power Double Weibul ##############################################
#####################################################################

dpweibul <- function(x, c,lambda){
  lambda*((c/2)*abs(x)^(c-1)*exp(-abs(x)^c))*F_dweibul(x, c)^(lambda-1)
}

c = 2
x <- seq(-10,10, by = 0.1)
lambda1 = dpweibul(x,c,lambda=1)
lambda2 = dpweibul(x,c,lambda=0.5)
lambda3 = dpweibul(x,c,lambda=4)

probs <- data.frame(x=x, Link = rep(c("\u03BB=1","\u03BB=0.5", "\u03BB=4"), each = length(x)),
                    probs = c(lambda1, lambda2, lambda3))
library(ggplot2)
ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDA")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F12_powerlomax.png", width = 11, height = 7)

#######################
### Distribuição Acumulada
#######################

Fpweibul <- function(x, c,lambda){
  F_dweibul(x, c)^lambda
}
c=2
### Gráfico
x <- seq(-10,10, by = 0.1)
lambda1 = Fpweibul(x,lambda=1, c=c)
lambda2 = Fpweibul(x,lambda=0.5, c=c)
lambda3 = Fpweibul(x,lambda=2, c=c)
lambda4 = Fpweibul(x,lambda=4, c=c)

probs <- data.frame(x=x, Link = rep(c("\u03BB=1","\u03BB=0.5", "\u03BB=2","\u03BB=4"), each = length(x)),
                    probs = c(lambda1, lambda2, lambda3, lambda4))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDA")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F13_Fpowerlomax.png", width = 11, height = 7)

##############################################3
##### Proporção de 0 e 1 ##################
#############################################
n = 1000
r = 100
beta0 = 0
beta1 = 1
lambda = c(0.5, 0.75, 1, 2, 4)

c= 2
prop = data.frame(matrix(ncol=2, nrow=0))
colnames(prop) = c("prop", "lambda")
for (j in 1:5){  
  for (i in 1:r){
    set.seed(i)
    x  = runif(1000, -10,10)
    eta = beta0 + beta1 * x
    p <- Fpweibul(eta, c,lambda[j])
    y = rbinom(1000,1,p)
    prop[i+r*(j-1),] = c(mean(y), lambda[j])
  }
}
prop$lambda = factor(prop$lambda)

require(ggplot2)
ggplot(prop, aes(y=prop, x=lambda)) + 
  geom_boxplot()+ylab('Proporções')+xlab('Potência')+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F17_propPlomax.png", width = 11, height = 7)



###################################################################################
##### REversa de potencia Lomax #############################
#################################################################
dinvpweibul <- function(x, lambda){
  lambda*((c/2)*abs(x)^(c-1)*exp(-abs(x)^c))*F_dweibul(-x, c)^(lambda-1)
}

x <- seq(-10,10, by = 0.1)
lambda1 = dinvpweibul(x,lambda=1)
lambda2 = dinvpweibul(x,lambda=0.5)
lambda3 = dinvpweibul(x,lambda=4)

probs <- data.frame(x=x, Link = rep(c("\u03BB=1","\u03BB=0.5", "\u03BB=4"), each = length(x)),
                    probs = c(lambda1, lambda2, lambda3))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDA")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F14_rpowerlomax.png", width = 11, height = 7)

#######################
### Distribuição Acumulada
#######################

Finvplomax <- function(x, lambda){
  1 - F_dlomax(-x)^lambda
}

## Simulação
r= 100
n=1000
parametros = data.frame(matrix(ncol=3,nrow=0))
colnames(parametros) = c("beta0", "beta1", "lambda")
beta0 = 0
beta1 = 1
lambda = 2

loglike_dpinvlomax <- function(x, y, param){
  eta <-x %*% param[-length(param)]
  p <- Finvplomax(eta, param[length(param)])
  logvero = log(prod(p^y * (1-p)^(1-y)))
  return(logvero) 
}


for (i in 1:r){
  set.seed(i)
  x  = runif(n, -10, 10)
  eta = beta0 + beta1 * x
  p <- Finvplomax(eta, lambda)
  y = rbinom(1000,1,p)
  dados <- data.frame("x" = x,"y" = y)
  
  m.log = try(glm(formula = y ~ x, family = binomial(link = "logit"), data=dados), TRUE)
  
  start = c(m.log$coefficients[1],m.log$coefficients[2], 1)
  estim = maxLik(logLik=loglike_dpinvlomax,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS")
  parametros[i,] = estim$estimate
}

colMeans(parametros)



### Gráfico
x <- seq(-10,10, by = 0.1)
lambda1 = Finvplomax(x,lambda=1)
lambda2 = Finvplomax(x,lambda=0.5)
lambda3 = Finvplomax(x,lambda=2)
lambda4 = Finvplomax(x,lambda=4)

probs <- data.frame(x=x, Link = rep(c("\u03BB=1","\u03BB=0.5", "\u03BB=2","\u03BB=4"), each = length(x)),
                    probs = c(lambda1, lambda2, lambda3, lambda4))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype = Link))+
  xlab("x")+ ylab("FDA")+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F15_Frpowerlomax.png", width = 11, height = 7)


##############################################3
##### Proporção de 0 e 1 ##################
#############################################
n = 1000
r = 100
beta0 = 0
beta1 = 1
lambda = c(1/5, 1/2, 1, 2, 3)

prop = data.frame(matrix(ncol=2, nrow=0))
colnames(prop) = c("prop", "lambda")
for (j in 1:5){  
  for (i in 1:r){
    set.seed(i)
    x  = runif(1000, -10,10)
    eta = beta0 + beta1 * x
    p <- Finvplomax(eta, lambda[j])
    y = rbinom(1000,1,p)
    prop[i+r*(j-1),] = c(mean(y), lambda[j])
  }
}
prop$lambda = factor(prop$lambda)

require(ggplot2)
ggplot(prop, aes(y=prop, x=lambda)) + 
  geom_boxplot()+ylab('Proporções')+xlab('Potência')+
  theme_bw()+ theme(legend.title = element_blank())
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F16_prop1_pLomax.png", width = 11, height = 7)




##################################33
###  Double Lomax ################
####################################
eta <- seq(-10,10, by = 0.1)

F_dlomax <- function(x){
  prob = vector(length = length(x))
  
  prob[which(x>0)] = 1 - 0.5*(1/(1+x[x>0]))
  prob[which(x<=0)] = 0.5*(1/(1 - x[x<=0]))
  
  return(prob)
}

plot(F_dlomax(eta))

r= 100
n=1000
parametros = data.frame(matrix(ncol=2,nrow=0))
colnames(parametros) = c("beta0", "beta1")
beta0 = 0
beta1 = 1

loglike_dlomax <- function(x, y, param){
  eta <-x %*% param
  p <- F_dlomax(eta)
  logvero = log(prod(p^y * (1-p)^(1-y)))
  return(logvero) 
}


for (i in 1:r){
  set.seed(i)
  x  = runif(n, -10, 10)
  eta = beta0 + beta1 * x
  p <- F_dlomax(eta)
  y = rbinom(1000,1,p)
  dados <- data.frame("x" = x,"y" = y)
  
  m.log = try(glm(formula = y ~ x, family = binomial(link = "logit"), data=dados), TRUE)
  
  start = c(m.log$coefficients[1],m.log$coefficients[2])
  estim = maxLik(logLik=loglike_dlomax,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS")
  parametros[i,] = estim$estimate
}

colMeans(parametros)

##########################################3
#### Plot lomax ################################33
#####################################################
library(ggplot2)

lomax <- function(x){1/(2*(1+abs(x))^2)}
laplace <- function(x){0.5* exp(-abs(x))}
### gráfico com as link functions
eta <- seq(-10,10, by = 0.1)

probs <- data.frame(eta = eta,`Distribuição` = rep(c("DLomax","Normal", "Cauchy", "Laplace"), each = length(eta)),
                    probs = c(lomax(eta), dnorm(eta), dcauchy(eta), laplace(eta)))

 ggplot(probs, aes(x=eta, y=probs, group = `Distribuição`)) +geom_line(aes(linetype=`Distribuição`, color = `Distribuição`))+
   scale_color_manual(values=c('black','red', 'black', 'black'))+
  xlab("x")+ ylab("FDP")+
  theme_bw()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "Lomax_comparacao.png", width = 6, height = 4)

### Comparação Link

### link functions
logit <- function(eta){exp(eta)/(1+exp(eta))}
probit <- function(eta){pnorm(eta)}

### gráfico com as link functions
eta <- seq(-15,15, by = 0.1)

probs <- data.frame(eta = eta, Link = rep(c("Logit","Probit", "DLomax"), each = length(eta)),
                    probs = c(logit(eta), probit(eta), F_dlomax(eta) ))

ggplot(probs, aes(x=eta, y=probs, group = Link)) +
  geom_line(aes(linetype=Link, color = Link))+
  scale_color_manual(values=c('red','black', 'black'))+
  xlab(expression(eta[i]))+ ylab(expression(p[i]))+
  theme_bw()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F9_lomax_link_comparation.png", width = 6, height = 4)


####################################################################
### Power Double Lomax ##############################################
#####################################################################

dplomax <- function(x, lambda){
  lambda*(0.5*1/(1+abs(x))^2)*F_dlomax(x)^(lambda-1)
}
plot(1/(1+abs(x))^2)

x <- seq(-10,10, by = 0.1)
lambda1 = dplomax(x,lambda=0.25)
lambda2 = dplomax(x,lambda=0.5)
lambda3 = dplomax(x,lambda=0.75)

plot(dplomax(x,lambda=1))

probs <- data.frame(x=x, Link = rep(c("\u03BB=0.25","\u03BB=0.5", "\u03BB=0.75"), each = length(x)),
                    probs = c(lambda1, lambda2, lambda3))
library(ggplot2)
ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype=Link, color = Link))+
 scale_color_manual(values=c('red','black', 'blue'))+
  xlab("x")+ ylab("FDP")+ylim(0,0.5)+
  theme_bw()+ theme(legend.title = element_blank(), legend.position = c(0.1, 0.85))
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "powerlomax_assimetriamenorque1.png", width = 6, height = 4)

#######################
### Distribuição Acumulada
#######################

Fplomax <- function(x, lambda){
  F_dlomax(x)^lambda
}
## Simulação
r= 100
n=1000
parametros = data.frame(matrix(ncol=3,nrow=0))
colnames(parametros) = c("beta0", "beta1", "lambda")
beta0 = 0
beta1 = 1
lambda = 2 

loglike_dplomax <- function(x, y, param){
  eta <-x %*% param[-length(param)]
  p <- Fplomax(eta, param[length(param)])
  logvero = log(prod(p^y * (1-p)^(1-y)))
  return(logvero) 
}


for (i in 1:r){
  set.seed(i)
  x  = runif(n, -10, 10)
  eta = beta0 + beta1 * x
  p <- Fplomax(eta, lambda)
  y = rbinom(1000,1,p)
  dados <- data.frame("x" = x,"y" = y)
  
  m.log = try(glm(formula = y ~ x, family = binomial(link = "logit"), data=dados), TRUE)
  
  start = c(m.log$coefficients[1],m.log$coefficients[2], 1)
  estim = maxLik(logLik=loglike_dplomax,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS")
  parametros[i,] = estim$estimate
}

colMeans(parametros)


### Gráfico
x <- seq(-10,10, by = 0.1)
lambda1 = Fplomax(x,lambda=0.25)
lambda2 = Fplomax(x,lambda=0.5)
lambda3 = Fplomax(x,lambda=1)
lambda4 = Fplomax(x,lambda=2)
lambda5 = Fplomax(x,lambda=6)

probs <- data.frame(x=x, Link = rep(c("\u03BB=0.25","\u03BB=0.5", "\u03BB=1","\u03BB=2", "\u03BB=6"), each = length(x)),
                    probs = c(lambda1, lambda2, lambda3, lambda4, lambda5))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype=Link, color = Link))+
  scale_color_manual(values=c('red', "orange", 'black',"green",'blue'))+
  xlab("x")+ ylab("FDA")+ ylim(0,1)+
  theme_bw()+ theme(legend.title = element_blank(), legend.position = c(0.1, 0.8))
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "powerlomaxacumulada_lambdas.png", width = 6, height = 4)


##############################################3
##### Proporção de 0 e 1 ##################
#############################################
n = 1000
r = 100
beta0 = 0
beta1 = 1
lambda = c(0.05,0.1,0.2,0.3,0.5, 1, 2, 4, 6, 8, 10, 12)

prop = data.frame(matrix(ncol=2, nrow=0))
colnames(prop) = c("prop", "lambda")
for (j in 1:11){  
  for (i in 1:r){
    set.seed(i)
    x  = runif(1000, -10,10)
    eta = beta0 + beta1 * x
    p <- Fplomax(eta, lambda[j])
    y = rbinom(1000,1,p)
    prop[i+r*(j-1),] = c(mean(y), lambda[j])
  }
}
prop$lambda = factor(prop$lambda)

require(ggplot2)
ggplot(prop, aes(y=prop, x=lambda)) + 
  geom_boxplot()+ylab('Proporções')+xlab('\u03BB')+
  theme_bw()+ theme(legend.title = element_blank())+ylim(0,1)
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F17_propPlomax.png", width = 6, height = 4)


mean(prop[prop$lambda==10,"prop"])
###################################################################################
##### REversa de potencia Lomax #############################
#################################################################
dinvplomax <- function(x, lambda){
  lambda*(0.5*1/(1+abs(x))^2)*F_dlomax(-x)^(lambda-1)
}

x <- seq(-10,10, by = 0.1)
lambda1 = dinvplomax(x,lambda=2)
lambda2 = dinvplomax(x,lambda=4)
lambda3 = dinvplomax(x,lambda=6)

plot(dplomax(x,lambda=1))

probs <- data.frame(x=x, Link = rep(c("\u03BB=2","\u03BB=4", "\u03BB=6"), each = length(x)),
                    probs = c(lambda1, lambda2, lambda3))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype=Link, color = Link))+
  scale_color_manual(values=c('red','black', 'blue'))+
  xlab("x")+ ylab("FDP")+ylim(0,0.5)+
  theme_bw()+ theme(legend.title = element_blank(), legend.position = c(0.9, 0.85))
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "rpowerlomaxmaiorque1.png", width = 6, height = 4)

#######################
### Distribuição Acumulada
#######################

Finvplomax <- function(x, lambda){
  1 - F_dlomax(-x)^lambda
}

## Simulação
r= 100
n=1000
parametros = data.frame(matrix(ncol=3,nrow=0))
colnames(parametros) = c("beta0", "beta1", "lambda")
beta0 = 0
beta1 = 1
lambda = 2

loglike_dpinvlomax <- function(x, y, param){
  eta <-x %*% param[-length(param)]
  p <- Finvplomax(eta, param[length(param)])
  logvero = log(prod(p^y * (1-p)^(1-y)))
  return(logvero) 
}


for (i in 1:r){
  set.seed(i)
  x  = runif(n, -10, 10)
  eta = beta0 + beta1 * x
  p <- Finvplomax(eta, lambda)
  y = rbinom(1000,1,p)
  dados <- data.frame("x" = x,"y" = y)
  
  m.log = try(glm(formula = y ~ x, family = binomial(link = "logit"), data=dados), TRUE)
  
  start = c(m.log$coefficients[1],m.log$coefficients[2], 1)
  estim = maxLik(logLik=loglike_dpinvlomax,start=start, y=dados[,2], x= model.matrix(~dados[,1]), method="BFGS")
  parametros[i,] = estim$estimate
}

colMeans(parametros)



### Gráfico
x <- seq(-10,10, by = 0.1)
lambda1 = Finvplomax(x,lambda=0.25)
lambda2 = Finvplomax(x,lambda=0.5)
lambda3 = Finvplomax(x,lambda=1)
lambda4 = Finvplomax(x,lambda=2)
lambda5 = Finvplomax(x,lambda=6)

probs <- data.frame(x=x, Link = rep(c("\u03BB=0.25","\u03BB=0.5", "\u03BB=1","\u03BB=2","\u03BB=6" ), each = length(x)),
                    probs = c(lambda1, lambda2, lambda3, lambda4, lambda5))

ggplot(probs, aes(x=x, y=probs, group = Link)) +geom_line(aes(linetype=Link, color = Link))+
  scale_color_manual(values=c('red', "orange", 'black',"green",'blue'))+
  xlab("x")+ ylab("FDA")+ ylim(0,1)+
  theme_bw()+ theme(legend.title = element_blank(), legend.position = c(0.1, 0.8))
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "Frpowerlomaxlambdas.png", width = 6, height = 4)


##############################################3
##### Proporção de 0 e 1 ##################
#############################################
n = 1000
r = 100
beta0 = 0
beta1 = 1
lambda = c(0.05,0.1,0.2,0.3,0.5, 1, 2, 4, 6, 8, 10, 12)

prop = data.frame(matrix(ncol=2, nrow=0))
colnames(prop) = c("prop", "lambda")
for (j in 1:11){  
  for (i in 1:r){
    set.seed(i)
    x  = runif(1000, -10,10)
    eta = beta0 + beta1 * x
    p <- Finvplomax(eta, lambda[j])
    y = rbinom(1000,1,p)
    prop[i+r*(j-1),] = c(mean(y), lambda[j])
  }
}
prop$lambda = factor(prop$lambda)

require(ggplot2)
ggplot(prop, aes(y=prop, x=lambda)) + 
  geom_boxplot()+ylab('Proporções')+xlab("\u03BB")+
  theme_bw()+ theme(legend.title = element_blank())+ylim(0,1)
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F16_prop1_pLomax.png", width = 6, height = 4)

mean(prop[prop$lambda==10,"prop"])