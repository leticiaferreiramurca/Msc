library(rstan)
options(mc.cores = 7)
library(loo)

source("funcoes_aplicacao.R")

####################################################################################
########### Blood Donation Dataset ##########################################################
###################################################################################
#https://archive.ics.uci.edu/ml/datasets/Blood+Transfusion+Service+Center

df = read.csv("/home/notmyname/Downloads/transfusion.data", header=T, fileEncoding = "utf8")
colnames(df) = sapply(strsplit(colnames(df),"[..]"), `[`, 1)

summary(df)

apply(df,2,max)
apply(df,2,sd) 

library(corrplot)

corrplot(cor(df[,-5]), is.corr = FALSE, method = 'color', col.lim = c(-1, 1),, col = colorRampPalette(c("blue","white","firebrick3"))(200), addCoef.col = 'grey30')

library(ggplot2)
tbl <- with(df, prop.table(table(whether)))
ggplot(as.data.frame(tbl), aes(factor(whether), Freq )) +
  scale_y_continuous(labels=scales::percent)+
geom_col(position = 'dodge')+ theme_bw()+ylab("Porcentagem")+xlab("")

summary(df)

Y = df[,5]
X = scale(df[,c(1,2,4)])


### Double Lomax
model_dlomax <- stan_model('dlomax.stan')
fit_dlomax <- sampling(model_dlomax, list(n = length(Y), y=Y, X = as.matrix(X), k = ncol(X)), iter = 5000, chains = 4, seed =1)
log_lik_1 <- extract_log_lik(fit_dlomax)
loo(log_lik_1)$looic 
waic(log_lik_1)$waic 

postSamples <- Reduce(cbind, extract(fit_dlomax, pars=c("beta0", "beta")))
dic(model.matrix(~X), Y, loglike_dlomax, postSamples)
eaic(model.matrix(~X), Y, loglike_dlomax, postSamples)
ebic(model.matrix(~X), Y, loglike_dlomax, postSamples)

### Reverse Power Lomax
model_rplomax <- stan_model('drplomax.stan')
fit_rplomax <- sampling(model_rplomax, list(n = length(Y), y=Y, X = as.matrix(X), k = ncol(X)), iter = 5000, chains = 4, seed =1)
log_lik_1 <- extract_log_lik(fit_rplomax)
loo(log_lik_1)$looic
waic(log_lik_1)$waic  

postSamples <- Reduce(cbind, extract(fit_rplomax, pars=c("beta0", "beta", "loglambda")))
dic(model.matrix(~X), Y, loglike_dpinvlomax, postSamples)
eaic(model.matrix(~X), Y, loglike_dpinvlomax, postSamples)
ebic(model.matrix(~X), Y, loglike_dpinvlomax, postSamples)

### Power Lomax
model_plomax <- stan_model('dplomax.stan')
fit_plomax <- sampling(model_plomax, list(n = length(Y), y=Y, X = as.matrix(X), k = ncol(X)), iter = 5000, chains = 4, seed = 1)
log_lik_1 <- extract_log_lik(fit_plomax)
loo(log_lik_1)$looic
waic(log_lik_1)$waic 

postSamples <- Reduce(cbind, extract(fit_plomax, pars=c("beta0", "beta", "loglambda")))
dic(model.matrix(~X), Y, loglike_dplomax, postSamples)
eaic(model.matrix(~X), Y, loglike_dplomax, postSamples)
ebic(model.matrix(~X), Y, loglike_dplomax, postSamples)

### Logistica 
model_logistic <- stan_model('dlogistic.stan')
fit_logistic <- sampling(model_logistic, list(n = length(Y), y=Y, X = as.matrix(X), k = ncol(X)), iter = 5000, chains = 4, seed =1)
log_lik_1 <- extract_log_lik(fit_logistic)
loo(log_lik_1)$looic
waic(log_lik_1)$waic 

postSamples <- Reduce(cbind, extract(fit_logistic, pars=c("beta0", "beta")))
dic(model.matrix(~X), Y, loglike_log, postSamples)
eaic(model.matrix(~X), Y, loglike_log, postSamples)
ebic(model.matrix(~X), Y, loglike_log, postSamples)

### Cloglog
model_cloglog <- stan_model('dcloglog.stan')
fit_cloglog <- sampling(model_cloglog, list(n = length(Y), y=Y, X = as.matrix(X), k = ncol(X)), iter = 5000, chains = 4, seed =1)
log_lik_1 <- extract_log_lik(fit_cloglog)
loo(log_lik_1)$looic
waic(log_lik_1)$waic 

postSamples <- Reduce(cbind, extract(fit_cloglog, pars=c("beta0", "beta")))
dic(model.matrix(~X), Y, loglike_cloglog, postSamples)
eaic(model.matrix(~X), Y, loglike_cloglog, postSamples)
ebic(model.matrix(~X), Y, loglike_cloglog, postSamples)

### Loglog
model_loglog <- stan_model('loglog.stan')
fit_loglog <- sampling(model_loglog, list(n = length(Y), y=Y, X = as.matrix(X), k = ncol(X)), iter = 5000, chains = 4, seed =1)
log_lik_1 <- extract_log_lik(fit_loglog)
loo(log_lik_1)$looic
waic(log_lik_1)$waic 

postSamples <- Reduce(cbind, extract(fit_loglog, pars=c("beta0", "beta")))
dic(model.matrix(~X), Y, loglike_loglog, postSamples)
eaic(model.matrix(~X), Y, loglike_loglog, postSamples)
ebic(model.matrix(~X), Y, loglike_loglog, postSamples)

### Probit
model_probit <- stan_model('probit.stan')
fit_probit <- sampling(model_probit, list(n = length(Y), y=Y, X = as.matrix(X), k = ncol(X)), iter = 5000, chains = 4, seed =1)
log_lik_1 <- extract_log_lik(fit_probit)
loo(log_lik_1)$looic
waic(log_lik_1)$waic

postSamples <- Reduce(cbind, extract(fit_probit, pars=c("beta0", "beta")))
dic(model.matrix(~X), Y, loglike_probit, postSamples)
eaic(model.matrix(~X), Y, loglike_probit, postSamples)
ebic(model.matrix(~X), Y, loglike_probit, postSamples)

### Cauchit
model_cauchit <- stan_model('cauchit.stan')
fit_cauchit <- sampling(model_cauchit, list(n = length(Y), y=Y, X = as.matrix(X), k = ncol(X)), iter = 5000, chains = 4, seed =1)
log_lik_1 <- extract_log_lik(fit_cauchit)
loo(log_lik_1)$looic
waic(log_lik_1)$waic 

postSamples <- Reduce(cbind, extract(fit_cauchit, pars=c("beta0", "beta")))
dic(model.matrix(~X), Y, loglike_cauchit, postSamples)
eaic(model.matrix(~X), Y, loglike_cauchit, postSamples)
ebic(model.matrix(~X), Y, loglike_cauchit, postSamples)

## tabela
modelos = list(fit_dlomax, fit_rplomax, fit_plomax, fit_logistic, fit_cloglog, fit_loglog, fit_probit,fit_cauchit)
loglikes = list(loglike_dlomax,  loglike_dpinvlomax, loglike_dplomax, loglike_log, loglike_cloglog, loglike_loglog,loglike_probit, loglike_cauchit )
medidas = data.frame(matrix(ncol = 9,nrow = length(modelos)))
colnames(medidas) = c("dist", "rho", "dbar", "dhat", "dic", "eaic", "ebic", "loo", "waic")
#modelos = list(fit_probit, fit_plomax)
#loglikes = list(loglike_probit,loglike_dplomax)
for(i in 1:length(modelos)){
  
  log_lik_1 <- extract_log_lik(modelos[[i]])
  medidas[i,"loo"] = loo(log_lik_1)$looic
  medidas[i,"waic"] = waic(log_lik_1)$waic 
  
  if(i %in% c(2,3)){
    postSamples <- Reduce(cbind, extract(modelos[[i]], pars=c("beta0", "beta", "loglambda")))
  }else{
    postSamples <- Reduce(cbind, extract(modelos[[i]], pars=c("beta0", "beta")))
  }
  
  dicc = dic(model.matrix(~X), Y, loglikes[[i]], postSamples)
  
  medidas[i,"dbar"] = dicc$dbar
  medidas[i,"dhat"] = dicc$dhat
  medidas[i,"rho"] = dicc$rho
  medidas[i,"dic"] = dicc$dic
  medidas[i,"eaic"] = eaic(model.matrix(~X), Y, loglikes[[i]], postSamples)
  medidas[i,"ebic"] = ebic(model.matrix(~X), Y, loglikes[[i]], postSamples)
}

medidas[,-1]
xtable(medidas[,-1], digits =3)

### valor dos parametros 

modelos = list(fit_dlomax, fit_rplomax, fit_plomax, fit_logistic, fit_cloglog, fit_loglog, fit_probit,fit_cauchit)
parametros = data.frame(matrix(ncol = 7,nrow = length(modelos)))
colnames(parametros) = c("dist", "beta0", "beta1", "beta2", "beta3", "loglambda", "lambda")


for (i in 1:length(modelos)){
  parametros[i, "beta0"] = mean(extract(modelos[[i]], permuted = TRUE)$beta0)
  parametros[i, "beta1"] = mean(extract(modelos[[i]], permuted = TRUE)$beta[,1])
  parametros[i, "beta2"] = mean(extract(modelos[[i]], permuted = TRUE)$beta[,2])
  parametros[i, "beta3"] = mean(extract(modelos[[i]], permuted = TRUE)$beta[,3])
  if(i %in% c(2,3)){
    parametros[i, "loglambda"] = mean(extract(modelos[[i]], permuted = TRUE)$loglambda)
    parametros[i, "lambda"] = mean(exp(extract(modelos[[i]], permuted = TRUE)$loglambda))
  }
  print(i)
}
library(xtable)
xtable(round(parametros,3), row.names = F, digits = 3)


rplomaxx <- extract(fit_rplomax,permuted = T)

par1 <- data.frame(matrix(ncol = 5))
names(par1) <- c("par", "mean", "sd", "median", "5", "95")

par1[1,1] <- mean(rplomaxx$beta0)
par1[1,2] <-sd(rplomaxx$beta0)
par1[1,3] <-median(rplomaxx$beta0)
par1[1,4:5] <- quantile(rplomaxx$beta0, c(0.05, 0.95))

par1[2:4,1] <-apply(rplomaxx$beta, 2, mean)
par1[2:4,2] <-apply(rplomaxx$beta, 2, sd)
par1[2:4,3] <-apply(rplomaxx$beta, 2, median)
par1[2:4,4:5] <-t(apply(rplomaxx$beta, 2, quantile, c(0.05, 0.95)))


par1[5,1] <-mean(exp(rplomaxx$loglambda))
par1[5,2] <- sd(exp(rplomaxx$loglambda))
par1[5,3] <-median(exp(rplomaxx$loglambda))
par1[5,4:5] <-quantile(exp(rplomaxx$loglambda), c(0.05, 0.95))
xtable(par1, digits = 3)

#########################################################################3
#### Análise Preditiva #####################################################
######################################################################
library(pROC)
library("cutpointr")
library(caret)

loglog <- function(x){exp(-exp(-x))}
cloglog <- function(x){1-exp(-exp(x))}

med_pred = data.frame(matrix(ncol = 7,nrow = length(modelos)))
names(med_pred) = c("cp","AUC", "Acuracia", "Precisao", "F1", "SENS", "ESP")
probabilidades <- list(F_dlomax, Finvplomax, Fplomax, plogis, cloglog, loglog, pnorm, pcauchy)

for (i in 1:nrow(parametros)){
  xbeta <- model.matrix(~X) %*% as.numeric(parametros[i,2:5],nrow=1)
  
  if(i %in% c(2,3)){
    prob <- sapply(xbeta, probabilidades[[i]],  loglambda = parametros[i,6])
  }else{
    prob <- sapply(xbeta, probabilidades[[i]])
  }
  
  cp <- 0.5
  
  med_pred[i,"cp"] <- cp
  prediction <- ifelse(prob>cp, 1,0)
  
  roc_obj <- pROC::roc(Y,prediction)
  
  med_pred[i,"AUC"] = pROC::auc(roc_obj)
  
  xtab <- table(prediction, Y)
  print(xtab)

  cm <- confusionMatrix(xtab, mode = "everything", positive = "1")
  if(i ==1){
    confusion <- xtab
  }else{
    confusion <- rbind(confusion,xtab) 
  }
  
  med_pred[i,"Acuracia"] <- cm$overall['Accuracy']
  med_pred[i, c("Precisao", "F1", "SENS", "ESP")] <- cm$byClass[c("Precision", "F1","Sensitivity","Specificity" )]
print(i)
}
library(xtable) 
xtable(round(med_pred,3), row.names = F, digits = 3)
xtable(confusion)
############################################################################
### Residuals ##########################################################
##########################################################################
## usei i=2
xbeta <- model.matrix(~X) %*% as.numeric(parametros[2,2:5],nrow=1)

prob <- Finvplomax(xbeta, parametros[2,6])

set.seed(1)
a <- pbinom(Y-1, 1, prob) 
b <- pbinom(Y, 1, prob)
res <- qnorm(runif(length(a),a, b))
plot(res)
qqnorm(res, pch = 1, frame = FALSE)
qqline(res, col = "steelblue", lwd = 2)
hist(res)

require(qqplotr)
require(qqplotr)
ggplot(data = data.frame(res), mapping = aes(sample = res)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Quantis Teóricos", y = "Resíduos Quantílicos")+theme_bw()
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "QQplot1.png", width = 6, height = 4)

gg <- ggplot(data = data.frame(res), aes(x = res)) +
  geom_histogram(colour="black", fill="white", bins =15) +
  labs(x = "Resíduos Quantílicos", y = "Frequência")+theme_bw()
gg
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "histres1.png", width = 6, height = 4)


data = data.frame("res" = res, "index" = 1:length(res))
gg <- ggplot(data = data, aes(x = index, y = res)) +
  geom_point(size = 2) +
  labs(x = "Index", y = "Resíduos Quantílicos")+theme_bw()
gg
ggsave(path = '/home/notmyname/Desktop/mestrado/imagens', filename = "quantres1.png", width = 6, height = 4)
