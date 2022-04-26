rm(list = ls(all = TRUE))

library(igraph)
library(spatialprobit)
library(corrplot)
source('/Users/mark/Documents/networks/scripts/gibbs_zsamplingfunctions.R')
source('/Users/mark/Documents/networks/scripts/sar_base.R')
source('/Users/mark/Documents/networks/scripts/sarprobit_sarcontinuous_functions.R')

setwd("C:\\Users\\mark\\Downloads\\ToreOpsahl")
#setwd("/data/blackgull/roeling/bayes")

W = as(as.matrix(readMM(file='Wmatrix.01.sym.rowstochastic.txt')), "dgCMatrix") 
#W = as(as.matrix(readMM(file='Wmatrix.meanmessages.rowstochastic.txt')), "dgCMatrix") 

FBattributes = read.table("attributes.csv", sep = ";", h = T)
FBattributes$gender01 = FBattributes$gender
FBattributes[FBattributes$gender01 == 1, ]$gender01 <- 0 # recode gender
FBattributes[FBattributes$gender01 == 2, ]$gender01 <- 1 # recode gender

FBattributes = FBattributes[order(FBattributes$id), ] # sort on ID to match W
# add intercept, remove ID and send.friends (equal to outdegree) 
FBattributes = as.data.frame(cbind(1, FBattributes[,-c(1, which(names(FBattributes) == "send.friends"))])) 
colnames(FBattributes)[1] = "Intercept"
head(FBattributes)

###### read in for testing popular #####
X = FBattributes[,c(1,10,3,8,9)]
y = FBattributes[,6] # popular
beta = list(); beta[[1]] = rep(0, ncol(X)) # for 1 y variable
#beta = rep(0, ncol(X))
computeMarginalEffects=TRUE
showProgress=TRUE

prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1)

prior2 = c()
prior2$novi = 1

#Sys.time()
#model1 = sar_continuous_mcmc(log(y), X, W, ndraw=100000, burn.in=1000, thinning=1, 
#                             prior=list(a1=1, a2=1, c=rep(0, ncol(X)), 
#                                    T=diag(ncol(X))*1e12, lflag = 0), 
#                             start=list(rho=0.75, beta=rep(0, ncol(X))), 
#                             m=10, computeMarginalEffects=TRUE, showProgress=TRUE)
#Sys.time()
#write(paste(model1$beta, model1$beta_std, model1$rho), "model1nomissing.txt")

#lmmodel1 = lm(log(popular) ~ gender01 + yearofstudy + averageChar + dayactive + send.day34th.friends, data = FBattributes)

# compare betas for log(popular)
#model1$beta
#summary(lmmodel1)

# add 1 
#W = as(as.matrix(W + 1), "dgCMatrix")

# simulate missing data
set.seed(19985)
missing.sample.10 = sample(nrow(FBattributes), 190)

set.seed(19985)
missing.sample.25 = sample(nrow(FBattributes), 475)

set.seed(19985)
missing.sample.35 = sample(nrow(FBattributes), 664)

set.seed(19985)
missing.sample.50 = sample(nrow(FBattributes), 950)

#yf1 = y; yf1[missing.sample.1] <- NA
yf10 = y; yf10[missing.sample.10] <- NA
yf25 = y; yf25[missing.sample.25] <- NA
yf35 = y; yf35[missing.sample.35] <- NA
yf50 = y; yf50[missing.sample.50] <- NA

#### cut model
impmodel.10 = sar_combined_mcmc_imp_cut(y = log(yf10), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"), model = 10)

impmodel.25 = sar_combined_mcmc_imp_cut(y = log(yf25), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"), model = 25)

impmodel.35 = sar_combined_mcmc_imp_cut(y = log(yf35), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"), model = 35)

Sys.time()
impmodel.50 = sar_combined_mcmc_imp_cut(y = log(yf50), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"), model = 50)
Sys.time()

y = FBattributes[,6] # popular

means.values.1 = colMeans(impmodel.1$y_imp[[1]])
means.values.10 = colMeans(impmodel.10$y_imp[[1]])
means.values.25 = colMeans(impmodel.25$y_imp[[1]])
means.values.35 = colMeans(impmodel.35$y_imp[[1]])

error.estimates.1 = as.data.frame(cbind(log(y[missing.sample.1]), means.values.1))
error.estimates.10 = as.data.frame(cbind(log(y[missing.sample.10]), means.values.10))
error.estimates.25 = as.data.frame(cbind(log(y[missing.sample.25]), means.values.25))
error.estimates.35 = as.data.frame(cbind(log(y[missing.sample.35]), means.values.35))

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

write(paste(impmodel.1$beta, impmodel.1$beta_std, impmodel.1$rho, mean(error.estimates.1$error), sd(error.estimates.1$error)), "cutmodel1.txt")
write(paste(impmodel.10$beta, impmodel.10$beta_std, impmodel.10$rho, mean(error.estimates.10$error), sd(error.estimates.10$error)), "cutmodel10.txt")
write(paste(impmodel.25$beta, impmodel.25$beta_std, impmodel.25$rho, mean(error.estimates.25error), sd(error.estimates.25$error)), "cutmodel25.txt")
write(paste(impmodel.35$beta, impmodel.35$beta_std, impmodel.35$rho, mean(error.estimates.35$error), sd(error.estimates.35$error)), "cutmodel35.txt")

#### full bayes 
impmodel.1fb = sar_combined_mcmc_imp_fbs(y = log(yf1), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                       prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                       m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"), model = 10)

#### full bayes 
impmodel.10fb = sar_combined_mcmc_imp_fbs(y = log(yf10), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                         prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                         m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"), model = 10)

#### full bayes 
impmodel.25fb = sar_combined_mcmc_imp_fbs(y = log(yf25), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                         prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                         m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"), model = 25)

impmodel.35fb = sar_combined_mcmc_imp_fbs(y = log(yf35), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                          prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                          m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"), model = 35)

Sys.time()
impmodel.50fb = sar_combined_mcmc_imp_fbs(y = log(yf50), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                          prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                          m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"), model = 50)
Sys.time()

y = FBattributes[,6] # popular

means.values.1fb = colMeans(impmodel.1fb$y_imp[[1]])
means.values.10fb = colMeans(impmodel.10fb$y_imp[[1]])
means.values.25fb = colMeans(impmodel.25fb$y_imp[[1]])
means.values.35fb = colMeans(impmodel.35fb$y_imp[[1]])
error.estimates.1 = as.data.frame(cbind(log(y[missing.sample.1]), means.values.1fb))
error.estimates.10 = as.data.frame(cbind(log(y[missing.sample.10]), means.values.10fb))
error.estimates.25 = as.data.frame(cbind(log(y[missing.sample.25]), means.values.25fb))
error.estimates.35 = as.data.frame(cbind(log(y[missing.sample.35]), means.values.35fb))
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

write(paste(impmodel.1fb$beta, impmodel.1fb$beta_std, impmodel.1fb$rho, mean(error.estimates.1$error), sd(error.estimates.1$error)), "fbsmodel1.txt")
write(paste(impmodel.10fb$beta, impmodel.10fb$beta_std, impmodel.10fb$rho, mean(error.estimates.10$error), sd(error.estimates.10$error)), "fbsmodel10.txt")
write(paste(impmodel.25fb$beta, impmodel.25fb$beta_std, impmodel.25fb$rho, mean(error.estimates.25$error), sd(error.estimates.25$error)), "fbsmodel25.txt")
write(paste(impmodel.35fb$beta, impmodel.35fb$beta_std, impmodel.35fb$rho, mean(error.estimates.35$error), sd(error.estimates.35$error)), "fbsmodel35.txt")

###### read in for testing gender #####
X = FBattributes[,c(1,3,5,6,7,8,9)]
# X = FBattributes[,c(1,3,5)]
y = FBattributes[,10] # gender
beta = list(); beta[[1]] = rep(0, ncol(X)) # for 1 y variable
#beta = rep(0, ncol(X))
computeMarginalEffects=TRUE
showProgress=TRUE

prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1)

prior2 = c()
prior2$novi = 1

# runned this function in the function script manually
model2 = sarprobit(y, X, W, ndraw=10000, burn.in=1000, thinning=1, 
                   prior=list(a1=1, a2=1, c=rep(0, ncol(X)), 
                   T=diag(ncol(X))*1e12, lflag = 0), 
                   start=list(rho=0.39, beta=rep(0, ncol(X))), 
                   m=10, computeMarginalEffects=TRUE, showProgress=TRUE)

lmmodel2 = lm(log(popular) ~ gender + yearofstudy + averageChar + dayactive + send.day34th.friends, data = FBattributes)



rhoestimate = read.table("/Users/mark/Documents/networks/rhoestimate50cut.txt")
rhoestimate = rhoestimate[c(1:100000),]
plot(rhoestimate, type = "l", xlab = "Rho", ylab = "iterations", main = "cut 50% random missingess")

rhoestimate = read.table("/Users/mark/Documents/networks/rhoestimate50fbs.txt")
rhoestimate = rhoestimate[c(1:100000),]
plot(rhoestimate, type = "l", xlab = "Rho", ylab = "iterations", main = "fbs 50 random missingess")



