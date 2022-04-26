rm(list = ls(all = TRUE))

library(igraph)
library(spatialprobit)
library(corrplot)
library(truncnorm)
library(snowboot)

source('D:/Documents/networks/scripts/functions/gibbs_zsamplingfunctions.R')
source('D:/Documents/networks/scripts/functions/sar_base.R')
source('D:/Documents/networks/scripts/functions/sar_continuous.R')

setwd("D:/Documents/downloads/ToreOpsahl")

OClinks_w = read.table("OClinks_w.txt")
colnames(OClinks_w) = c("id1", "id2", "messages")
OClinks_w_g = graph_from_edgelist(as.matrix(OClinks_w[,c(1:2)]), directed = F)
net = igraph_to_network(OClinks_w_g)

#set.seed(19985)
#missing.sample.1 = LSMI(net, n.seeds = 5, n.neigh = 1, seeds = NULL, classic = F)
#length(missing.sample.1$sampleN)
set.seed(19986)
missing.sample.10 = LSMI(net, n.seeds = 2, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.10$sampleN))/1899

set.seed(19986)
missing.sample.25 = LSMI(net, n.seeds = 46, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.25$sampleN))/1899

set.seed(19986)
missing.sample.35 = LSMI(net, n.seeds = 80, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.35$sampleN))/1899

set.seed(19986)
missing.sample.50 = LSMI(net, n.seeds = 184, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.50$sampleN))/1899

W = as(as.matrix(readMM(file='Wmatrix.01.sym.rowstochastic.txt')), "dgCMatrix") 

FBattributes = read.table("attributes.csv", sep = ";", h = T)
FBattributes$gender01 = FBattributes$gender
FBattributes[FBattributes$gender01 == 1, ]$gender01 <- 0 # recode gender
FBattributes[FBattributes$gender01 == 2, ]$gender01 <- 1 # recode gender

FBattributes = FBattributes[order(FBattributes$id), ] # sort on ID to match W
# add intercept, remove ID and send.friends (equal to outdegree) 
FBattributes = as.data.frame(cbind(1, FBattributes[,-c(1, which(names(FBattributes) == "send.friends"))])) 
colnames(FBattributes)[1] = "Intercept"
head(FBattributes)

###### read in for testing Year of Study #####
X = as.matrix(FBattributes[,c(1,10,8,9)])
y = FBattributes[,3] # year of study 
#y[y > 3] <- 3
y = as.integer(y)

model1ord = sar_ordered_probit_mcmc(y = y, X = X, W = W, showProgress=TRUE)


#################################################
#################################################

beta = list(); beta[[1]] = rep(0, ncol(X)) # for 1 y variable
#beta = rep(0, ncol(X))
computeMarginalEffects=TRUE
showProgress=TRUE

prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1)
start=list(rho=0.39, rhosd=0.16, beta=rep(0, ncol(X)), phi=c(-Inf, 0:(max(y)-1), Inf))

prior2 = c()
prior2$novi = 1

#yf1 = y; yf1[unique(sort(missing.sample.1$sampleN))] <- NA
yf10 = y; yf10[unique(sort(missing.sample.10$sampleN))] <- NA
yf25 = y; yf25[unique(sort(missing.sample.25$sampleN))] <- NA
yf35 = y; yf35[unique(sort(missing.sample.35$sampleN))] <- NA
yf50 = y; yf50[unique(sort(missing.sample.50$sampleN))] <- NA

#### cut model
source("D:/Documents/networks/scripts/functions/sar_imp_cut.R")
impmodel.10 = sar_combined_mcmc_imp_cut(y = yf10, x=X, W, ndraw=1000, burn.in=10, thinning=1,start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, model = 10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("orderedprobit"))

impmodel.25 = sar_combined_mcmc_imp_cut(y = yf25, x=X, W, ndraw=100000, burn.in=1000, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, model = 25, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("orderedprobit"))

Sys.time()
impmodel.35 = sar_combined_mcmc_imp_cut(y = yf35, x=X, W, ndraw=100000, burn.in=1000, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, model = 35, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("orderedprobit"))
Sys.time()


#y[5] <- NA
model1 = sar_combined_mcmc_imp_fbs(y = yf1, x = X, W = W,
                                   method = c("orderedprobit"),
                                   thinning=1,
                                   prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                   m=10,
                                   computeMarginalEffects=TRUE,
                                   showProgress=TRUE,
                                   start = start)

model10 = sar_combined_mcmc_imp_fbs(y = yf10, x = X, W = W,
                                   method = c("orderedprobit"),
                                   thinning=1,
                                   prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                   m=10,
                                   computeMarginalEffects=TRUE,
                                   showProgress=TRUE,
                                   start = start)

model25 = sar_combined_mcmc_imp_fbs(y = yf25, x = X, W = W,
                                   method = c("orderedprobit"),
                                   thinning=1,
                                   prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                   m=10,
                                   computeMarginalEffects=TRUE,
                                   showProgress=TRUE,
                                   start = start)

model35 = sar_combined_mcmc_imp_fbs(y = yf35, x = X, W = W,
                                    method = c("orderedprobit"),
                                    thinning=1,
                                    prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                    m=10,
                                    computeMarginalEffects=TRUE,
                                    showProgress=TRUE,
                                    start = start)

y = FBattributes[,3] # YOS

means.values.1fb = colMeans(impmodel.1fb$y_imp[[1]])
means.values.10fb = colMeans(impmodel.10fb$y_imp[[1]])
means.values.25fb = colMeans(impmodel.25fb$y_imp[[1]])
means.values.35fb = colMeans(impmodel.35fb$y_imp[[1]])
error.estimates.1 = as.data.frame(cbind(log(y[unique(sort(missing.sample.1$sampleN))]), means.values.1fb))
error.estimates.10 = as.data.frame(cbind(log(y[unique(sort(missing.sample.10$sampleN))]), means.values.10fb))
error.estimates.25 = as.data.frame(cbind(log(y[unique(sort(missing.sample.25$sampleN))]), means.values.25fb))
error.estimates.35 = as.data.frame(cbind(log(y[unique(sort(missing.sample.35$sampleN))]), means.values.35fb))
error.estimates.1$error = abs(error.estimates.1[,1] - error.estimates.1[,2])
error.estimates.10$error = abs(error.estimates.10[,1] - error.estimates.10[,2])
error.estimates.25$error = abs(error.estimates.25[,1] - error.estimates.25[,2])
error.estimates.35$error = abs(error.estimates.35[,1] - error.estimates.35[,2])
mean(error.estimates.1$error)
mean(error.estimates.10$error)
mean(error.estimates.25$error)
mean(error.estimates.35$error)
sd(error.estimates.1$error)
sd(error.estimates.10$error)
sd(error.estimates.25$error)
sd(error.estimates.35$error)
write(paste(impmodel.1fb$beta, impmodel.1fb$beta_std, impmodel.1fb$rho, mean(error.estimates.1$error), sd(error.estimates.1$error)), "fbsmodel1YOSsnow.txt")
write(paste(impmodel.10fb$beta, impmodel.10fb$beta_std, impmodel.10fb$rho, mean(error.estimates.10$error), sd(error.estimates.10$error)), "fbsmodel10YOSsnow.txt")
write(paste(impmodel.25fb$beta, impmodel.25fb$beta_std, impmodel.25fb$rho, mean(error.estimates.25$error), sd(error.estimates.25$error)), "fbsmodel25YOSsnow.txt")
write(paste(impmodel.35fb$beta, impmodel.35fb$beta_std, impmodel.35fb$rho, mean(error.estimates.35$error), sd(error.estimates.35$error)), "fbsmodel35YOSsnow.txt")

imputed1 = impmodel.1fb$y_imp[[1]]
imputed10 = impmodel.10fb$y_imp[[1]]
imputed25 = impmodel.25fb$y_imp[[1]]
imputed35 = impmodel.35fb$y_imp[[1]]
write.table(imputed1, "imputed_cutmodel_YOS1.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed10, "imputed_cutmodel_YOS10.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed25, "imputed_cutmodel_YOS25.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed35, "imputed_cutmodel_YOS35.csv", sep = ";", quote = F, col.names=F, row.names=F)

beta1 = impmodel.1fb$b
beta10 = impmodel.10fb$b
beta25 = impmodel.25fb$b
beta35 = impmodel.35fb$b

write.table(beta1, "beta_cutmodel_YOS1.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(beta10, "beta_cutmodel_YOS10.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(beta25, "beta_cutmodel_YOS25.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(beta35, "beta_cutmodel_YOS35.csv", sep = ";", quote = F, col.names=F, row.names=F)


