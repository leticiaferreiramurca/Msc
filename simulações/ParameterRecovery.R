################################################################################
###### Parameter Recovery ######################################################
################################################################################
library(rstan)
options(mc.cores = 8)

vies_comparacao <- data.frame(matrix(nrow = 0, ncol = 6))
colnames(vies_comparacao) <-  c("modelo","n","beta0", "beta1", "lambda", "lambdai")
rmse_comparacao <- data.frame(matrix(nrow = 0, ncol = 6))
colnames(rmse_comparacao) <-  c("modelo","n","beta0", "beta1", "lambda","lambdai")

########################################################
#### Estimação Logistica ###########################
#########################################################
#Generate Fake data
model_logistic <- stan_model('dlogistic.stan')
n = c(500, 1000, 2000)
r = 100 
k=0
beta0 = 0 
beta1 = 1

estimation <-data.frame(matrix(ncol=7, nrow=r*length(n)))
colnames(estimation) <- c("n", "beta0", "upperBeta0", "lowerBeta0", "beta1", "upperBeta1", "lowerBeta1")

start_time <- Sys.time()
### 100 amostras
for(ni in n){  
  XX = matrix(ncol = r, nrow = ni)
  yy = matrix(ncol = r, nrow = ni)
  for(i in 1:r){
    set.seed(i)
    XX[,i] = runif(ni, -3, 3)
    eta = beta0 + beta1*XX[,i]
    p <- plogis(eta) 
    yy[,i] = rbinom(ni,1,p)
  }
  
  #Estimation
  
  for(i in 1:r){
    fit_lomax <- sampling(model_logistic, list(n = ni, y=yy[,i], x=XX[,i]), iter = 1000, chains = 4)
    
    estimation[k+i,1] <- ni
    estimation[k+i,2] <- mean(extract(fit_lomax, permuted = TRUE)$beta0)
    estimation[k+i,5] <- mean(extract(fit_lomax, permuted = TRUE)$beta1)
    
    estimation[k+i,3] <- as.data.frame(summary(fit_lomax)[[1]])["beta0","2.5%"]
    estimation[k+i,4] <- as.data.frame(summary(fit_lomax)[[1]])["beta0","97.5%"]
    
    estimation[k+i, 6] <- as.data.frame(summary(fit_lomax)[[1]])["beta1","2.5%"]
    estimation[k+i, 7] <- as.data.frame(summary(fit_lomax)[[1]])["beta1","97.5%"]
    print(k+i)
  }
  k = k + 100
}
end_time <- Sys.time()

end_time - start_time

#write.csv(estimation, "estimation_logistic.csv")

logistic <- read.csv( "estimation_logistic.csv")[,-1]

beta0 = 0
beta1 = 1
vies <- data.frame(matrix(ncol = 3))
names(vies) <- c("n", "beta0", "beta1")
rmse <- data.frame(matrix(ncol = 3))
names(rmse) <- c("n", "beta0", "beta1")
pc <- vector()
alpha1 <- vector()
alpha2 <- vector()
n = c(500, 1000, 2000)

## Viés
#beta0
i=1
for(ni in n){
  print( ni)
  estimation_ni = estimation[estimation$n == ni,]
  ##Viés
  print("viés")
  vies[i,2] <- mean(estimation_ni[,2]) - beta0
  vies[i,3] <-  mean(estimation_ni[,5]) - beta1
  vies[i,1] <- ni
  
  ##rmse
  print("rmse")
  rmse[i,2] <- mean((estimation_ni[,2] - beta0)^2)
  rmse[i,3] <- mean((estimation_ni[,5] - beta1)^2)
  rmse[i,1] <- ni
  
  ## probabilidade de cobertura
  print("PC")
  print(mean(beta0 > estimation_ni[,3] & beta0 < estimation_ni[,4]))
  print(mean(beta1 > estimation_ni[,6] & beta1 < estimation_ni[,7]))
  
  #alpha1
  print("alpha1")
  print(mean(beta0 < estimation_ni[,3]))
  print(mean(beta1 < estimation_ni[,6]))
  
  #alpha2
  print("alpha2")
  print(mean(beta0 > estimation_ni[,4]))
  print(mean(beta1 > estimation_ni[,7]))
  i = i+1
}
nrow1 = nrow(vies)
vies_comparacao[1:nrow1,2:4] = vies
vies_comparacao[1:nrow1,1] = rep("Logistica", nrow1)

rmse_comparacao[1:nrow1,2:4] = rmse
rmse_comparacao[1:nrow1,1] = rep("Logistica", nrow1)
########################################################
#### Estimação Double Lomax ###########################
#########################################################
#Generate Fake data
F_dlomax <- function(x){
  prob = vector(length = length(x))
  
  prob[which(x>0)] = 1 - 0.5*(1/(1+x[x>0]))
  prob[which(x<=0)] = 0.5*(1/(1 - x[x<=0]))
  
  return(prob)
}

model_lomax <- stan_model('dlomax.stan')
n = c(500, 1000, 2000)
r = 100 
k=0

estimation <-data.frame(matrix(ncol=7, nrow=r*length(n)))
colnames(estimation) <- c("n", "beta0", "upperBeta0", "lowerBeta0", "beta1", "upperBeta1", "lowerBeta1")

start_time <- Sys.time()
### 100 amostras
for(ni in n){  
  XX = matrix(ncol = r, nrow = ni)
  yy = matrix(ncol = r, nrow = ni)
  for(i in 1:r){
    set.seed(i)
    XX[,i] = runif(ni, -3, 3)
    eta = beta0 + beta1*XX[,i]
    p <- F_dlomax(eta) 
    yy[,i] = rbinom(ni,1,p)
  }
  
  #Estimation
  
  for(i in 1:r){
    fit_lomax <- sampling(model_lomax, list(n = ni, y=yy[,i], x=XX[,i]), iter = 200, chains = 4)
    
    estimation[k+i,1] <- ni
    estimation[k+i,2] <- mean(extract(fit_lomax, permuted = TRUE)$beta0)
    estimation[k+i,5] <- mean(extract(fit_lomax, permuted = TRUE)$beta1)
    
    estimation[k+i,3] <- as.data.frame(summary(fit_lomax)[[1]])["beta0","2.5%"]
    estimation[k+i,4] <- as.data.frame(summary(fit_lomax)[[1]])["beta0","97.5%"]
    
    estimation[k+i, 6] <- as.data.frame(summary(fit_lomax)[[1]])["beta1","2.5%"]
    estimation[k+i, 7] <- as.data.frame(summary(fit_lomax)[[1]])["beta1","97.5%"]
    print(k+i)
  }
  k = k + 100
}
end_time <- Sys.time()

end_time - start_time

dlomax <- read.csv( "estimation_dlomax.csv")[,-1]


vies <- data.frame(matrix(ncol = 3))
names(vies) <- c("n", "beta0", "beta1")
rmse <- data.frame(matrix(ncol = 3))
names(rmse) <- c("n", "beta0", "beta1")
pc <- vector()
alpha1 <- vector()
alpha2 <- vector()
i=1
## Viés
#beta0
for(ni in n){
  print( ni)
  estimation_ni = estimation[estimation$n == ni,]
  ##Viés
  ##Viés
  print("viés")
  vies[i,2] <- mean(estimation_ni[,2]) - beta0
  vies[i,3] <-  mean(estimation_ni[,5]) - beta1
  vies[i,1] <- ni
  
  ##rmse
  print("rmse")
  rmse[i,2] <- mean((estimation_ni[,2] - beta0)^2)
  rmse[i,3] <- mean((estimation_ni[,5] - beta1)^2)
  rmse[i,1] <- ni
  
  ## probabilidade de cobertura
  print("PC")
  print(mean(beta0 > estimation_ni[,3] & beta0 < estimation_ni[,4]))
  print(mean(beta1 > estimation_ni[,6] & beta1 < estimation_ni[,7]))
  
  #alpha1
  print("alpha1")
  print(mean(beta0 < estimation_ni[,3]))
  print(mean(beta1 < estimation_ni[,6]))
  
  #alpha2
  print("alpha2")
  print(mean(beta0 > estimation_ni[,4]))
  print(mean(beta1 > estimation_ni[,7]))
  i=i+1
}

nrow2 = nrow(vies)+nrow1
vies_comparacao[(nrow1+1):nrow2,2:4] = vies
vies_comparacao[(nrow1+1):nrow2,1] = rep("DLomax", nrow(vies))

rmse_comparacao[(nrow1+1):nrow2,2:4] = rmse
rmse_comparacao[(nrow1+1):nrow2,1] = rep("DLomax", nrow(rmse))

##############################################
### Estimação Potência Lomax ####################################################
###################################################################################
#generate fake data
Fplomax <- function(x, lambda){
  F_dlomax(x)^lambda
}
model_plomax <- stan_model('dplomax.stan')

lambdas <- c(0.25,0.5, 2, 4)
n = c(500, 1000, 2000)
r = 100 
k=0

estimation <-data.frame(matrix(ncol=11, nrow=r*length(n)*length(lambdas)))
colnames(estimation) <- c("n","lambda", "beta0", "upperBeta0", "lowerBeta0", "beta1", "upperBeta1", "lowerBeta1", "lambda","upperlambda", "lowerlambda")


for(lambdai in lambdas){
  for(ni in n){  
    XX = matrix(ncol = r, nrow = ni)
    yy = matrix(ncol = r, nrow = ni)
    for(i in 1:r){
      set.seed(i)
      XX[,i] = runif(ni, -3, 3)
      eta = beta0 + beta1*XX[,i]
      p <- Fplomax(eta, lambdai) 
      yy[,i] = rbinom(ni,1,p)
    }
    
    #Estimation
    
    for(i in 1:r){
      fit_plomax <- sampling(model_plomax, list(n = ni, y=yy[,i], x=XX[,i]), iter = 200, chains = 4)
      
      estimation[k+i,1] <- ni
      estimation[k+i,2] <- lambdai 
      estimation[k+i,3] <- mean(extract(fit_plomax, permuted = TRUE)$beta0)
      estimation[k+i,6] <- mean(extract(fit_plomax, permuted = TRUE)$beta1)
      estimation[k+i,9] <- mean(extract(fit_plomax, permuted = TRUE)$loglambda)
      
      estimation[k+i,4] <- as.data.frame(summary(fit_plomax)[[1]])["beta0","2.5%"]
      estimation[k+i,5] <- as.data.frame(summary(fit_plomax)[[1]])["beta0","97.5%"]
      
      estimation[k+i, 7] <- as.data.frame(summary(fit_plomax)[[1]])["beta1","2.5%"]
      estimation[k+i, 8] <- as.data.frame(summary(fit_plomax)[[1]])["beta1","97.5%"]
      
      estimation[k+i, 10] <- as.data.frame(summary(fit_plomax)[[1]])["loglambda","2.5%"]
      estimation[k+i, 11] <- as.data.frame(summary(fit_plomax)[[1]])["loglambda","97.5%"]
      
      print(k+i)
    }
    k = k+100
  }
}
#write.csv(estimation, "estimation_plomax.csv")

plomax <- read.csv( "estimation_plomax.csv")[,-1]


lambdas <- c(0.25,0.5, 2, 4)
beta0 = 0 
beta1 = 1
vies <- data.frame(matrix(ncol = 5))

names(vies) <- c("n", "lambdai","beta0", "beta1", "lambda")
rmse <- data.frame(matrix(ncol = 5))
names(rmse) <- c("n", "lambdai", "beta0", "beta1", "lambda")
pc <- vector()
alpha1 <- vector()
alpha2 <- vector()
i=1

## Viés
#beta0
for (lambdai in lambdas){
  print("-------------------------------------------------------------")
  print(lambdai)
  print("-------------------------------------------------------------")
  for(ni in n){
    print( ni)
    estimation_ni = estimation[((estimation$n == ni) & (estimation$lambda == lambdai)),]
    
    ##Viés
    print("viés")
    vies[i,3] <- mean(estimation_ni[,3]) - beta0
    vies[i,4] <- mean(estimation_ni[,6]) - beta1
    vies[i,5] <- exp(mean(estimation_ni[,9])) - lambdai
    vies[i,1] <- ni
    vies[i,2] <- lambdai
    
    ##rmse
    print("rmse")
    rmse[i,3] <- mean((estimation_ni[,3] - beta0)^2)
    rmse[i,4] <- mean((estimation_ni[,6] - beta1)^2)
    rmse[i,5] <- mean((exp(estimation_ni[,9]) - lambdai)^2)
    rmse[i,1] <- ni
    rmse[i,2] <- lambdai
    
    
    ## probabilidade de cobertura
    print("PC")
    print(mean(beta0 > estimation_ni[,4] & beta0 < estimation_ni[,5]))
    print(mean(beta1 > estimation_ni[,7] & beta1 < estimation_ni[,8]))
    print(mean(lambdai > exp(estimation_ni[,10]) & lambdai < exp(estimation_ni[,11])))
    
    #alpha1
    print("alpha1")
    print(mean(beta0 < estimation_ni[,4]))
    print(mean(beta1 < estimation_ni[,7]))
    print(mean(lambdai < exp(estimation_ni[,10])))
    
    #alpha2
    print("alpha2")
    print(mean(beta0 > estimation_ni[,5]))
    print(mean(beta1 > estimation_ni[,8]))
    print(mean(lambdai > exp(estimation_ni[,11])))
    i = i+1
  }
}

nrow3 = nrow(vies)+nrow2
vies_comparacao[(nrow2+1):nrow3,2:5] = vies[,-2]
vies_comparacao[(nrow2+1):nrow3,1] = paste("PDL(", vies$lambdai,")", sep = "")
vies_comparacao[(nrow2+1):nrow3,6] = vies[,2]

rmse_comparacao[(nrow2+1):nrow3,2:5] = rmse[,-2]
rmse_comparacao[(nrow2+1):nrow3,1] = paste("PDL(", rmse$lambdai,")", sep = "")
rmse_comparacao[(nrow2+1):nrow3,6] = rmse[,2]

##############################################################################
### Estimação Reversa de Potência Lomax ####################################################
###################################################################################
#Generate Fake data
Finvplomax <- function(x, lambda){
  1 - F_dlomax(-x)^lambda
}

model_rplomax <- stan_model('drplomax.stan')

lambdas <- c(0.25,0.5, 2, 4)
n = c(500, 1000, 2000)
r = 100 
k=0

estimation <-data.frame(matrix(ncol=11, nrow=r*length(n)*length(lambdas)))
colnames(estimation) <- c("n","lambda", "beta0", "upperBeta0", "lowerBeta0", "beta1", "upperBeta1", "lowerBeta1", "lambda","upperlambda", "lowerlambda")


for(lambdai in lambdas){
  for(ni in n){  
    XX = matrix(ncol = r, nrow = ni)
    yy = matrix(ncol = r, nrow = ni)
    for(i in 1:r){
      set.seed(i)
      XX[,i] = runif(ni, -3, 3)
      eta = beta0 + beta1*XX[,i]
      p <- Finvplomax(eta, lambdai) 
      yy[,i] = rbinom(ni,1,p)
    }
    
    #Estimation
    
    for(i in 1:r){
      fit_rplomax <- sampling(model_rplomax, list(n = ni, y=yy[,i], x=XX[,i]), iter = 200, chains = 4)
      
      estimation[k+i,1] <- ni
      estimation[k+i,2] <- lambdai 
      estimation[k+i,3] <- mean(extract(fit_rplomax, permuted = TRUE) $beta0)
      estimation[k+i,6] <- mean(extract(fit_rplomax, permuted = TRUE)$beta1)  
      estimation[k+i,9] <- mean(extract(fit_rplomax, permuted = TRUE)$loglambda)
      
      estimation[k+i,4] <- as.data.frame(summary(fit_rplomax)[[1]])["beta0","2.5%"]
      estimation[k+i,5] <- as.data.frame(summary(fit_rplomax)[[1]])["beta0","97.5%"]
      
      estimation[k+i, 7] <- as.data.frame(summary(fit_rplomax)[[1]])["beta1","2.5%"]
      estimation[k+i, 8] <- as.data.frame(summary(fit_rplomax)[[1]])["beta1","97.5%"]
      
      estimation[k+i, 10] <- as.data.frame(summary(fit_rplomax)[[1]])["loglambda","2.5%"]
      estimation[k+i, 11] <- as.data.frame(summary(fit_rplomax)[[1]])["loglambda","97.5%"]
      
      print(k+i)
    }
    k = k+100
  }
}

#write.csv(estimation, "estimation_rplomax.csv")


rplomax <- read.csv( "estimation_rplomax.csv")[,-1]


lambdas <- c(0.25,0.5, 2, 4)
beta0 = 0 
beta1 = 1
vies <- data.frame(matrix(ncol = 5))

names(vies) <- c("n", "lambdai","beta0", "beta1", "lambda")
rmse <- data.frame(matrix(ncol = 5))
names(rmse) <- c("n", "lambdai", "beta0", "beta1", "lambda")
pc <- vector()
alpha1 <- vector()
alpha2 <- vector()
i=1


## Viés
#beta0
for (lambdai in lambdas){
  print("-------------------------------------------------------------")
  print(lambdai)
  print("-------------------------------------------------------------")
  for(ni in n){
    print( ni)
    estimation_ni = estimation[((estimation$n == ni) & (estimation$lambda == lambdai)),]
    ##Viés
    print("viés")
    vies[i,3] <- mean(estimation_ni[,3]) - beta0
    vies[i,4] <- mean(estimation_ni[,6]) - beta1
    vies[i,5] <- exp(mean(estimation_ni[,9])) - lambdai
    vies[i,1] <- ni
    vies[i,2] <- lambdai
    
    ##rmse
    print("rmse")
    rmse[i,3] <- mean((estimation_ni[,3] - beta0)^2)
    rmse[i,4] <- mean((estimation_ni[,6] - beta1)^2)
    rmse[i,5] <- mean((exp(estimation_ni[,9]) - lambdai)^2)
    rmse[i,1] <- ni
    rmse[i,2] <- lambdai
    
    ## probabilidade de cobertura
    print("PC")
    print(mean(beta0 > estimation_ni[,4] & beta0 < estimation_ni[,5]))
    print(mean(beta1 > estimation_ni[,7] & beta1 < estimation_ni[,8]))
    print(mean(lambdai > exp(estimation_ni[,10]) & lambdai < exp(estimation_ni[,11])))
    
    #alpha1
    print("alpha1")
    print(mean(beta0 < estimation_ni[,4]))
    print(mean(beta1 < estimation_ni[,7]))
    print(mean(lambdai < exp(estimation_ni[,10])))
    
    #alpha2
    print("alpha2")
    print(mean(beta0 > estimation_ni[,5]))
    print(mean(beta1 > estimation_ni[,8]))
    print(mean(lambdai > exp(estimation_ni[,11])))
    i=i+1
  }
}

nrow4 = nrow(vies)+nrow3
vies_comparacao[(nrow3+1):nrow4,2:5] = vies[,-2]
vies_comparacao[(nrow3+1):nrow4,1] = paste("RPDL(", vies$lambdai,")", sep = "")
vies_comparacao[(nrow3+1):nrow4,6] = vies[,2]

rmse_comparacao[(nrow3+1):nrow4,2:5] = rmse[,-2]
rmse_comparacao[(nrow3+1):nrow4,1] = paste("RPDL(", rmse$lambdai,")", sep = "")
rmse_comparacao[(nrow3+1):nrow4,6] = rmse[,2]

rmse_comparacao

#########################################################################
####### GRÁFICOS ######################################################
########################################################################
library("ggsci")
######################
## Viés
## beta 0 
ggplot(vies_comparacao, aes(x=n, y=beta0, group = modelo)) +geom_line(aes(color=modelo))+
  geom_point(aes(color=modelo))+
  ylab(expression(paste("Viés(",beta[0],")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()

ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F30_biasBeta0.png", width = 11, height = 7)

## beta 1
ggplot(vies_comparacao, aes(x=n, y=beta1, group = modelo)) +geom_line(aes(color=modelo))+
  geom_point(aes(color=modelo))+
  ylab(expression(paste("Viés(",beta[1],")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()

ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F31_biasBeta1.png", width = 11, height = 7)


## lambda
ggplot(na.omit(vies_comparacao), aes(x=n, y=lambda, group = modelo)) +geom_line(aes(color=modelo))+
  geom_point(aes(color=modelo))+
  ylab(expression(paste("Viés(",lambda,")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()

ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F32_biasLambda.png", width = 11, height = 7)


## EQM
## beta 0 
ggplot(rmse_comparacao, aes(x=n, y=beta0, group = modelo)) +geom_line(aes(color=modelo))+
  geom_point(aes(color=modelo))+
  ylab(expression(paste("EQM(",beta[0],")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F33_eqmBeta0.png", width = 11, height = 7)


## beta 1
ggplot(rmse_comparacao, aes(x=n, y=beta1, group = modelo)) +geom_line(aes(color=modelo))+
  geom_point(aes(color=modelo))+
  ylab(expression(paste("EQM(",beta[1],")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F34_eqmBeta1.png", width = 11, height = 7)


## lambda
ggplot(na.omit(rmse_comparacao), aes(x=n, y=lambda, group = modelo)) +geom_line(aes(color=modelo))+
  geom_point(aes(color=modelo))+
  ylab(expression(paste("RMSE(",lambda,")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F35_eqmLambda.png", width = 11, height = 7)

##########################################################
######################## plots ###########################
#########################################################
library(dplyr)

dfs <- list(logistic, dlomax, plomax, rplomax)
namesdist <- c("Logística", "DLomax", "PDLomax(", "RPDLomax(")

vies <- function(x, par){ mean(x) - par}
rms <- function(x, par){mean((x - par)^2)}
df_vieses <- data.frame(matrix(ncol = 5))
names(df_vieses) = c("dist", "n", "par","lambda", "vies")
df_rmses <- data.frame(matrix(ncol = 5))
names(df_rmses) = c("dist", "n", "par", "lambda", "rms")
k=1
for(i in 1:4){
  if(i <3){
    df_vieses[k:(k+2),c("n", "vies")] <- dfs[[i]] %>% group_by(n) %>% do(data.frame(val=vies(.$beta0, 0)))
    df_vieses[k:(k+2),"par"] <- "beta0"
    df_vieses[k:(k+2),"dist"] <- namesdist[i]
    k=k+3
    df_vieses[k:(k+2),c("n", "vies")] <- dfs[[i]] %>% group_by(n) %>% do(data.frame(val=vies(.$beta1, 1)))
    df_vieses[k:(k+2),"par"] <- "beta1"
    df_vieses[k:(k+2),"dist"] <- namesdist[i]
    k=k+3
  } else{
    df_vieses[k:(k+11),c("n", "lambda","vies")] <- dfs[[i]] %>% group_by(n,lambda) %>% do(data.frame(val=vies(.$beta0, 0)))
    df_vieses[k:(k+11),"par"] <- "beta0"
    df_vieses[k:(k+11),"dist"] <- paste(namesdist[i], df_vieses[k:(k+11),c("lambda")], ")", sep = "") 
    k=k+12
    df_vieses[k:(k+11),c("n", "lambda","vies")] <- dfs[[i]] %>% group_by(n,lambda) %>% do(data.frame(val=vies(.$beta1, 1)))
    df_vieses[k:(k+11),"par"] <- "beta1"
    df_vieses[k:(k+11),"dist"] <- paste(namesdist[i], df_vieses[k:(k+11),c("lambda")], ")", sep = "") 
    k=k+12
    df_vieses[k:(k+11),c("n", "lambda","vies")] <- dfs[[i]] %>% group_by(n,lambda) %>% do(data.frame(val=vies(exp(.$lambda.1),unique(.$lambda))))
    df_vieses[k:(k+11),"par"] <- "lambda"
    df_vieses[k:(k+11),"dist"] <- paste(namesdist[i], df_vieses[k:(k+11),c("lambda")], ")", sep = "") 
    k=k+12
  }
}

## rmses
for(i in 1:4){
  if(i <3){
    df_rmses[k:(k+2),c("n", "rms")] <- dfs[[i]] %>% group_by(n) %>% do(data.frame(val=rms(.$beta0, 0)))
    df_rmses[k:(k+2),"par"] <- "beta0"
    df_rmses[k:(k+2),"dist"] <- namesdist[i]
    k=k+3
    df_rmses[k:(k+2),c("n", "rms")] <- dfs[[i]] %>% group_by(n) %>% do(data.frame(val=rms(.$beta1, 1)))
    df_rmses[k:(k+2),"par"] <- "beta1"
    df_rmses[k:(k+2),"dist"] <- namesdist[i]
    k=k+3
  } else{
    df_rmses[k:(k+11),c("n", "lambda","rms")] <- dfs[[i]] %>% group_by(n,lambda) %>% do(data.frame(val=rms(.$beta0, 0)))
    df_rmses[k:(k+11),"par"] <- "beta0"
    df_rmses[k:(k+11),"dist"] <- paste(namesdist[i], df_rmses[k:(k+11),c("lambda")], ")", sep = "") 
    k=k+12
    df_rmses[k:(k+11),c("n", "lambda","rms")] <- dfs[[i]] %>% group_by(n,lambda) %>% do(data.frame(val=rms(.$beta1, 1)))
    df_rmses[k:(k+11),"par"] <- "beta1"
    df_rmses[k:(k+11),"dist"] <- paste(namesdist[i], df_rmses[k:(k+11),c("lambda")], ")", sep = "") 
    k=k+12
    df_rmses[k:(k+11),c("n", "lambda","rms")] <- dfs[[i]] %>% group_by(n,lambda) %>% do(data.frame(val=rms(exp(.$lambda.1),unique(.$lambda))))
    df_rmses[k:(k+11),"par"] <- "lambda"
    df_rmses[k:(k+11),"dist"] <- paste(namesdist[i], df_rmses[k:(k+11),c("lambda")], ")", sep = "") 
    k=k+12
  }
}
## Viés
library(ggplot2)
## beta 0 
vies_comparacao = df_vieses[df_vieses$par == "beta0",]
ggplot(vies_comparacao, aes(x=n, y=vies, group = dist)) +geom_line(aes(color=dist))+
  geom_point(aes(color=dist))+
  ylab(expression(paste("Viés(",beta[0],")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()

ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F30_biasBeta0.png", width = 6, height = 4)

## beta 1
vies_comparacao = df_vieses[df_vieses$par == "beta1",]
ggplot(vies_comparacao, aes(x=n, y=vies, group = dist)) +geom_line(aes(color=dist))+
  geom_point(aes(color=dist))+
  ylab(expression(paste("Viés(",beta[1],")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()

ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F31_biasBeta1.png", width = 6, height = 4)


## lambda

vies_comparacao = df_vieses[df_vieses$par == "lambda",]
ggplot(vies_comparacao, aes(x=n, y=vies, group = dist)) +geom_line(aes(color=dist))+
  geom_point(aes(color=dist))+
  ylab(expression(paste("Viés(",lambda,")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()

ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F32_biasLambda.png", width = 6, height = 4)


## EQM
## beta 0 
rmse_comparacao = df_rmses[df_rmses$par == "beta0",]
ggplot(rmse_comparacao, aes(x=n, y=rms, group = dist)) +geom_line(aes(color=dist))+
  geom_point(aes(color=dist))+
  ylab(expression(paste("EQM(",beta[0],")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F33_eqmBeta0.png", width = 6, height = 4)


## beta 1
rmse_comparacao = df_rmses[df_rmses$par == "beta1",]
ggplot(rmse_comparacao, aes(x=n, y=rms, group = dist)) +geom_line(aes(color=dist))+
  geom_point(aes(color=dist))+
  ylab(expression(paste("EQM(",beta[1],")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F34_eqmBeta1.png", width = 6, height = 4)


## lambda
rmse_comparacao = df_rmses[df_rmses$par == "lambda",]
ggplot(rmse_comparacao, aes(x=n, y=rms, group = dist)) +geom_line(aes(color=dist))+
  geom_point(aes(color=dist))+
  ylab(expression(paste("EQM(",lambda,")")))+
  labs(color = "Link")+
  theme_bw()+scale_color_jco()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "F35_eqmLambda.png", width = 6, height = 4)
