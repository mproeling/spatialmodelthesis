rm(list = ls(all = TRUE))

library(igraph)
library(spatialprobit)
library(corrplot)
library(snowboot)
source('D:/Documents/networks/scripts/functions/gibbs_zsamplingfunctions.R')
source('D:/Documents/networks/scripts/functions/sar_base.R')
source('D:/Documents/networks/scripts/functions/sar_continuous.R')

setwd("D:/Documents/downloads/ToreOpsahl")
#setwd("/data/blackgull/roeling/bayes")

# snowball sampling to decide where to create the missing data
# example
tr <- make_tree(40, children = 3, mode = "undirected")
plot.igraph(tr)
net = igraph_to_network(tr)
LSMI(net, n.seeds = 3, n.neigh = 1, seeds = NULL, classic = F)

## 
OClinks_w = read.table("OClinks_w.txt")
colnames(OClinks_w) = c("id1", "id2", "messages")
OClinks_w_g = graph_from_edgelist(as.matrix(OClinks_w[,c(1:2)]), directed = F)
net = igraph_to_network(OClinks_w_g)

set.seed(19985)
missing.sample.10 = LSMI(net, n.seeds = 10, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.10$sampleN))

set.seed(19985)
missing.sample.25 = LSMI(net, n.seeds = 58, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.25$sampleN))

set.seed(19985)
missing.sample.35 = LSMI(net, n.seeds = 73, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.35$sampleN))/1899

set.seed(19985)
missing.sample.50 = LSMI(net, n.seeds = 159, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.50$sampleN))/1899
 
set.seed(19985)
missing.sample.75 = LSMI(net, n.seeds = 457, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.75$sampleN))/1899


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

Sys.time()
# this is an adapted version (for continuous y) of the sarprobit function from the spatialprobit package
model1 = sar_continuous_mcmc(log(y), X, W, ndraw=100000, burn.in=1000, thinning=1, 
                             prior=list(a1=1, a2=1, c=rep(0, ncol(X)), 
                                        T=diag(ncol(X))*1e12, lflag = 0), 
                             start=list(rho=0.75, beta=rep(0, ncol(X))), 
                             m=10, computeMarginalEffects=TRUE, showProgress=TRUE)
Sys.time()
write(paste(model1$beta, model1$beta_std, model1$rho), "model1nomissing.txt")

lmmodel1 = lm(log(popular) ~ gender01 + yearofstudy + averageChar + dayactive + send.day34th.friends, data = FBattributes)

# compare betas for log(popular)
model1$beta
summary(lmmodel1)

# add 1 
#W = as(as.matrix(W + 1), "dgCMatrix")

# simulate missing data
#yf1 = y; yf1[unique(sort(missing.sample.1$sampleN))] <- NA
yf10 = y; yf10[unique(sort(missing.sample.10$sampleN))] <- NA
yf25 = y; yf25[unique(sort(missing.sample.25$sampleN))] <- NA
yf35 = y; yf35[unique(sort(missing.sample.35$sampleN))] <- NA
yf50 = y; yf50[unique(sort(missing.sample.50$sampleN))] <- NA
yf75 = y; yf75[unique(sort(missing.sample.75$sampleN))] <- NA

#y[6]<-NA
#yf1 = y 
#### cut model
source("D:/Documents/networks/scripts/functions/sar_imp_cut.R")
impmodel.10 = sar_combined_mcmc_imp_cut(y = log(yf10), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                        model = 10)

impmodel.25 = sar_combined_mcmc_imp_cut(y = log(yf25), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                        model = 25)

impmodel.35 = sar_combined_mcmc_imp_cut(y = log(yf35), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                        model = 35)

impmodel.50 = sar_combined_mcmc_imp_cut(y = log(yf50), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                        model = 50)

impmodel.50 = sar_combined_mcmc_imp_cut(y = log(yf75), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                        model = 75)

#rhoall = read.table("rhoestimate_cut.txt")
#rhoall = rhoall[11001:22000,]

y = FBattributes[,6] # popular

write.table(dim(impmodel.1$y_imp[[1]]), "yimp.model1.csv", sep = ";", quote = F, row.names = F, col.names = F)
write.table(dim(impmodel.10$y_imp[[1]]), "yimp.model10.csv", sep = ";", quote = F, row.names = F, col.names = F)
write.table(dim(impmodel.25$y_imp[[1]]), "yimp.model25.csv", sep = ";", quote = F, row.names = F, col.names = F)
write.table(dim(impmodel.35$y_imp[[1]]), "yimp.model35.csv", sep = ";", quote = F, row.names = F, col.names = F)
write.table(dim(impmodel.50$y_imp[[1]]), "yimp.model50.csv", sep = ";", quote = F, row.names = F, col.names = F)

means.values.1 = colMeans(impmodel.1$y_imp[[1]])
means.values.10 = colMeans(impmodel.10$y_imp[[1]])
means.values.25 = colMeans(impmodel.25$y_imp[[1]])
means.values.35 = colMeans(impmodel.35$y_imp[[1]])
means.values.50 = colMeans(impmodel.50$y_imp[[1]])

error.estimates.1 = as.data.frame(cbind(log(y[unique(sort(missing.sample.1$sampleN))]), means.values.1))
error.estimates.10 = as.data.frame(cbind(log(y[unique(sort(missing.sample.10$sampleN))]), means.values.10))
error.estimates.25 = as.data.frame(cbind(log(y[unique(sort(missing.sample.25$sampleN))]), means.values.25))
error.estimates.35 = as.data.frame(cbind(log(y[unique(sort(missing.sample.35$sampleN))]), means.values.35))
error.estimates.50 = as.data.frame(cbind(log(y[unique(sort(missing.sample.50$sampleN))]), means.values.50))

error.estimates.1$error = abs(error.estimates.1[,1] - error.estimates.1[,2])
error.estimates.10$error = abs(error.estimates.10[,1] - error.estimates.10[,2])
error.estimates.25$error = abs(error.estimates.25[,1] - error.estimates.25[,2])
error.estimates.35$error = abs(error.estimates.35[,1] - error.estimates.35[,2])
error.estimates.50$error = abs(error.estimates.50[,1] - error.estimates.50[,2])

mean(error.estimates.1$error)
mean(error.estimates.10$error)
mean(error.estimates.25$error)
mean(error.estimates.35$error)
mean(error.estimates.50$error)
sd(error.estimates.1$error)
sd(error.estimates.10$error)
sd(error.estimates.25$error)
sd(error.estimates.35$error)
sd(error.estimates.50$error)

write(paste(impmodel.1$beta, impmodel.1$beta_std, impmodel.1$rho, mean(error.estimates.1$error), sd(error.estimates.1$error)), "cutmodel1snow.txt")
write(paste(impmodel.10$beta, impmodel.10$beta_std, impmodel.10$rho, mean(error.estimates.10$error), sd(error.estimates.10$error)), "cutmodel10snow.txt")
write(paste(impmodel.25$beta, impmodel.25$beta_std, impmodel.25$rho, mean(error.estimates.25$error), sd(error.estimates.25$error)), "cutmodel25snow.txt")
write(paste(impmodel.35$beta, impmodel.35$beta_std, impmodel.35$rho, mean(error.estimates.35$error), sd(error.estimates.35$error)), "cutmodel35snow.txt")
write(paste(impmodel.50$beta, impmodel.50$beta_std, impmodel.50$rho, mean(error.estimates.50$error), sd(error.estimates.50$error)), "cutmodel50snow.txt")

#### full bayes 
source("D:/Documents/networks/scripts/functions/sar_imp_fbs.R")
impmodel.1fb = sar_combined_mcmc_imp_fbs(y = log(yf1), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                       prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                       m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                       model = 1)

#### full bayes 
impmodel.10fb = sar_combined_mcmc_imp_fbs(y = log(yf10), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                         prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                         m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                         model = 10)

#### full bayes 
impmodel.25fb = sar_combined_mcmc_imp_fbs(y = log(yf25), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                         prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                         m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                         model = 25)

impmodel.35fb = sar_combined_mcmc_imp_fbs(y = log(yf35), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                          prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                          m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                          model = 35)

impmodel.50fb = sar_combined_mcmc_imp_fbs(y = log(yf50), x=X, W, ndraw=100000, burn.in=1000, thinning=1, 
                                          prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                          m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("continuous"),
                                          model = 50)

y = FBattributes[,6] # popular

write.table(dim(impmodel.1fb$y_imp[[1]]), "yimp.fbsmodel1.csv", sep = ";", quote = F, row.names = F, col.names = F)
write.table(dim(impmodel.10fb$y_imp[[1]]), "yimp.fbsmodel10.csv", sep = ";", quote = F, row.names = F, col.names = F)
write.table(dim(impmodel.25fb$y_imp[[1]]), "yimp.cutmodel25.csv", sep = ";", quote = F, row.names = F, col.names = F)
write.table(dim(impmodel.35fb$y_imp[[1]]), "yimp.cutmodel35.csv", sep = ";", quote = F, row.names = F, col.names = F)
write.table(dim(impmodel.50fb$y_imp[[1]]), "yimp.fbsmodel50.csv", sep = ";", quote = F, row.names = F, col.names = F)

means.values.1fb = colMeans(impmodel.1fb$y_imp[[1]])
means.values.10fb = colMeans(impmodel.10fb$y_imp[[1]])
means.values.25fb = colMeans(impmodel.25fb$y_imp[[1]])
means.values.35fb = colMeans(impmodel.35$y_imp[[1]])
means.values.50fb = colMeans(impmodel.50fb$y_imp[[1]])

error.estimates.1 = as.data.frame(cbind(log(y[unique(sort(missing.sample.1$sampleN))]), means.values.1fb))
error.estimates.10 = as.data.frame(cbind(log(y[unique(sort(missing.sample.10$sampleN))]), means.values.10fb))
error.estimates.25 = as.data.frame(cbind(log(y[unique(sort(missing.sample.25$sampleN))]), means.values.25fb))
error.estimates.35 = as.data.frame(cbind(log(y[unique(sort(missing.sample.35$sampleN))]), means.values.35fb))
error.estimates.50 = as.data.frame(cbind(log(y[unique(sort(missing.sample.50$sampleN))]), means.values.50fb))

error.estimates.1$error = abs(error.estimates.1[,1] - error.estimates.1[,2])
error.estimates.10$error = abs(error.estimates.10[,1] - error.estimates.10[,2])
error.estimates.25$error = abs(error.estimates.25[,1] - error.estimates.25[,2])
error.estimates.35$error = abs(error.estimates.35[,1] - error.estimates.35[,2])
error.estimates.50$error = abs(error.estimates.50[,1] - error.estimates.50[,2])

mean(error.estimates.1$error)
mean(error.estimates.10$error)
mean(error.estimates.25$error)
mean(error.estimates.35$error)
mean(error.estimates.50$error)

sd(error.estimates.1$error)
sd(error.estimates.10$error)
sd(error.estimates.25$error)
sd(error.estimates.35$error)
sd(error.estimates.50$error)

write(paste(impmodel.1fb$beta, impmodel.1fb$beta_std, impmodel.1fb$rho, mean(error.estimates.1$error), sd(error.estimates.1$error)), "fbsmodel1snow.txt")
write(paste(impmodel.10fb$beta, impmodel.10fb$beta_std, impmodel.10fb$rho, mean(error.estimates.10$error), sd(error.estimates.10$error)), "fbsmodel10snow.txt")
write(paste(impmodel.25fb$beta, impmodel.25fb$beta_std, impmodel.25fb$rho, mean(error.estimates.25$error), sd(error.estimates.25$error)), "fbsmodel25snow.txt")
write(paste(impmodel.35$beta, impmodel.35$beta_std, impmodel.35$rho, mean(error.estimates.35$error), sd(error.estimates.35$error)), "cutmodel35snow.txt")
write(paste(impmodel.50fb$beta, impmodel.50fb$beta_std, impmodel.50fb$rho, mean(error.estimates.50$error), sd(error.estimates.50$error)), "fbsmodel50snow.txt")

####################### 
# COMPARE # 
values = log(y[unique(missing.sample.1$sampleN)])
diff.cut = matrix(NA, nrow(impmodel.1$y_imp[[1]]), ncol(impmodel.1$y_imp[[1]]))
diff.fbs = matrix(NA, nrow(impmodel.1fb$y_imp[[1]]), ncol(impmodel.1fb$y_imp[[1]]))
for(i in 1:nrow(impmodel.1fb$y_imp[[1]])){
  diff.cut[i,] = values - impmodel.1$y_imp[[1]][i,]
  diff.fbs[i,] = values - impmodel.1fb$y_imp[[1]][i,]
}


impmodel.1fb$y_imp[[1]]
yf1 = y; yf1[unique(sort(missing.sample.1$sampleN))] <- NA
yf1[unique(sort(missing.sample.1$sampleN))] <- NA

# for the snowball sampled missing data
residuals = as.data.frame(cbind(log(y[unique(sort(missing.sample.10$sampleN))]), colMeans(impmodel.10$y_imp[[1]])))
residuals = as.data.frame(cbind(log(y[unique(sort(missing.sample.25$sampleN))]), colMeans(impmodel.25$y_imp[[1]])))
residuals = as.data.frame(cbind(log(y[unique(sort(missing.sample.35$sampleN))]), colMeans(impmodel.35$y_imp[[1]])))
residuals = as.data.frame(cbind(log(y[unique(sort(missing.sample.50$sampleN))]), colMeans(impmodel.50$y_imp[[1]])))

residuals = as.data.frame(cbind(log(y[unique(sort(missing.sample.10$sampleN))]), colMeans(impmodel.10fb$y_imp[[1]])))
residuals = as.data.frame(cbind(log(y[unique(sort(missing.sample.25$sampleN))]), colMeans(impmodel.25fb$y_imp[[1]])))
residuals = as.data.frame(cbind(log(y[unique(sort(missing.sample.35$sampleN))]), colMeans(impmodel.35fb$y_imp[[1]])))
residuals = as.data.frame(cbind(log(y[unique(sort(missing.sample.50$sampleN))]), colMeans(impmodel.50fb$y_imp[[1]])))

residuals$res = residuals$V1 - residuals$V2
residuals$res2 = residuals$res^2
mean(residuals$res)
mean(residuals$res^2)

write.table(residuals, "residuals.impmodel10", sep = ";", quote = F, col.names = T, row.names=F)
write.table(residuals, "residuals.impmodel25", sep = ";", quote = F, col.names = T, row.names=F)
write.table(residuals, "residuals.impmodel35", sep = ";", quote = F, col.names = T, row.names=F)
write.table(residuals, "residuals.impmodel50", sep = ";", quote = F, col.names = T, row.names=F)

write.table(impmodel.10$bdraw[[1]], "beta_cut10.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.25$bdraw[[1]], "beta_cut25.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.35$bdraw[[1]], "beta_cut35.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.50$bdraw[[1]], "beta_cut50.csv", sep = ";", quote = F, col.names = F, row.names = F)

write.table(impmodel.10$y_imp[[1]], "imp_cut10.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.25$y_imp[[1]], "imp_cut25.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.35$y_imp[[1]], "imp_cut35.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.50$y_imp[[1]], "imp_cut50.csv", sep = ";", quote = F, col.names = F, row.names = F)

save(objects,file="impmodel10cut_popularity.RData")

write.table(residuals, "residuals.impmodel10fb", sep = ";", quote = F, col.names = T, row.names=F)
write.table(residuals, "residuals.impmodel25fb", sep = ";", quote = F, col.names = T, row.names=F)
write.table(residuals, "residuals.impmodel35fb", sep = ";", quote = F, col.names = T, row.names=F)
write.table(residuals, "residuals.impmodel50fb", sep = ";", quote = F, col.names = T, row.names=F)

write.table(impmodel.10fb$bdraw[[1]], "beta_fbs10.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.25fb$bdraw[[1]], "beta_fbs25.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.35fb$bdraw[[1]], "beta_fbs35.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.50fb$bdraw[[1]], "beta_fbs50.csv", sep = ";", quote = F, col.names = F, row.names = F)

write.table(impmodel.10fb$y_imp[[1]], "imp_fbs10.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.25fb$y_imp[[1]], "imp_fbs25.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.35fb$y_imp[[1]], "imp_fbs35.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.50fb$y_imp[[1]], "imp_fbs50.csv", sep = ";", quote = F, col.names = F, row.names = F)

save(objects,file="impmodel50fbs_popularity.RData")






# for the completely random missing data:
residuals = as.data.frame(cbind(log(y[missing.sample.10]), colMeans(impmodel.10$y_imp[[1]])))
residuals = as.data.frame(cbind(log(y[missing.sample.25]), colMeans(impmodel.25$y_imp[[1]])))
residuals = as.data.frame(cbind(log(y[missing.sample.35]), colMeans(impmodel.35$y_imp[[1]])))
residuals = as.data.frame(cbind(log(y[missing.sample.50]), colMeans(impmodel.50$y_imp[[1]])))

