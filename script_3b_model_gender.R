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

set.seed(19985)
missing.sample.10 = LSMI(net, n.seeds = 10, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.10$sampleN))

set.seed(19985)
missing.sample.25 = LSMI(net, n.seeds = 58, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.25$sampleN))

set.seed(19985)
missing.sample.35 = LSMI(net, n.seeds = 73, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.35$sampleN))

set.seed(19985)
missing.sample.50 = LSMI(net, n.seeds = 159, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.50$sampleN))

set.seed(19985)
missing.sample.75 = LSMI(net, n.seeds = 457, n.neigh = 1, seeds = NULL, classic = F)
length(unique(missing.sample.75$sampleN))

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
X = as.matrix(FBattributes[,c(1,5,3,8)])
y = FBattributes[,10] # Gender 01 
#y[y > 3] <- 3
y = as.integer(y)

model1all = glm(y ~ as.matrix(FBattributes[,-c(1,2,10)]), family = "binomial")
model1gender = glm(y ~ X[,-1], family = "binomial")

#source("D:/Documents/networks/scripts/functions/sar_probit.R")
#model1gender = sar_probit_mcmc(y = y, X = X, W = W, ndraw=15000, burn.in=1000, showProgress=TRUE) # eigen versie met z zelf gemaakt
#model1gender = sar_probit_mcmc(y = y, X = as.matrix(FBattributes[,-c(2,6,10)]), W = W, ndraw=100, burn.in=10, showProgress=TRUE) # eigen versie met z zelf gemaakt


 
#model1gender = sar_combined_mcmc_imp_cut(y = y, x = X, W = W, showProgress=TRUE, method = "probit")
#model1gender = sar_combined_mcmc_imp_cut(y = y, x = X, W = W, ndraw=100000, burn.in=1000, showProgress=TRUE, method = "probit")
##################################################
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
yf75 = y; yf75[unique(sort(missing.sample.75$sampleN))] <- NA

#### cut model
source("D:/Documents/networks/scripts/functions/sar_imp_cut.R")
impmodel.10 = sar_combined_mcmc_imp_cut(y = yf10, x=X, W, ndraw=15000, burn.in=1000, thinning=1,start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = 10)

impmodel.25 = sar_combined_mcmc_imp_cut(y = yf25, x=X, W, ndraw=15000, burn.in=1000, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = 25)

Sys.time()
impmodel.35 = sar_combined_mcmc_imp_cut(y = yf35, x=X, W, ndraw=15000, burn.in=1000, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = 35)
Sys.time()

Sys.time()
impmodel.50 = sar_combined_mcmc_imp_cut(y = yf50, x=X, W, ndraw=15000, burn.in=1000, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = 50)
Sys.time()

impmodel.75 = sar_combined_mcmc_imp_cut(y = yf75, x=X, W, ndraw=15000, burn.in=1000, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = 75)


y = FBattributes[,10] # Gender01
imputed = impmodel.10$y_imp[[1]]
imputed = impmodel.25$y_imp[[1]]
imputed = impmodel.35$y_imp[[1]]
imputed = impmodel.50$y_imp[[1]]

write(paste(impmodel.1$beta, impmodel.1$beta_std, impmodel.1$rho, mean(error.estimates.1$error), sd(error.estimates.1$error)), "cutmodel1Gendersnow.txt")
write(paste(impmodel.10$beta, impmodel.10$beta_std, impmodel.10$rho, mean(error.estimates.10$error), sd(error.estimates.10$error)), "cutmodel10Gendersnow.txt")
write(paste(impmodel.25$beta, impmodel.25$beta_std, impmodel.25$rho, mean(error.estimates.25$error), sd(error.estimates.25$error)), "cutmodel25Gendersnow.txt")
write(paste(impmodel.35$beta, impmodel.35$beta_std, impmodel.35$rho, mean(error.estimates.35$error), sd(error.estimates.35$error)), "cutmodel35Gendersnow.txt")
write(paste(impmodel.50$beta, impmodel.50$beta_std, impmodel.50$rho, mean(error.estimates.50$error), sd(error.estimates.50$error)), "cutmodel50Gendersnow.txt")

save()

#y[5] <- NA
source("D:/Documents/networks/scripts/functions/sar_imp_fbs.R")
model10 = sar_combined_mcmc_imp_fbs(y = yf10, x = X, W = W, ndraw = 15000, burn.in = 1000, method = c("probit"),
                                   prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                   m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                   model = 10)

model25 = sar_combined_mcmc_imp_fbs(y = yf25, x = X, W = W, ndraw = 15000, burn.in = 1000, method = c("probit"),
                                   prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                   m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                   model = 25)

model35 = sar_combined_mcmc_imp_fbs(y = yf35, x = X, W = W, ndraw = 15000, burn.in = 1000, method = c("probit"),
                                    prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                    m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                    model = 35)

model50 = sar_combined_mcmc_imp_fbs(y = yf50, x = X, W = W, ndraw = 15000, burn.in = 1000, method = c("probit"),
                                    prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                    m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                    model = 50)

model75 = sar_combined_mcmc_imp_fbs(y = yf75, x = X, W = W, ndraw = 15000, burn.in = 1000, method = c("probit"),
                                    prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                    m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                    model = 75)


imputed10 = impmodel.10$y_imp[[1]]
imputed25 = impmodel.25$y_imp[[1]]
imputed35 = impmodel.35$y_imp[[1]]
imputed50 = impmodel.50$y_imp[[1]]
imputed75 = impmodel.75$y_imp[[1]]

write.table(imputed10, "imputed_cutmodel_gender10.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed25, "imputed_cutmodel_gender25.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed35, "imputed_cutmodel_gender35.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed50, "imputed_cutmodel_gender50.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed75, "imputed_cutmodel_gender75.csv", sep = ";", quote = F, col.names=F, row.names=F)

write.table(model1gender$bdraw, "bdraws_cutmodel_nomis.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.10$bdraw[[1]], "bdraws_cutmodel10.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.25$bdraw[[1]], "bdraws_cutmodel25.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.35$bdraw[[1]], "bdraws_cutmodel35.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.50$bdraw[[1]], "bdraws_cutmodel50.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.75$bdraw[[1]], "bdraws_cutmodel75.csv", sep = ";", quote = F, col.names = F, row.names = F)

save(objects,file="nomis_gender.RData")
save(objects,file="impmodel10cuts_gender.RData")
save(objects,file="impmodel25cuts_gender.RData")
save(objects,file="impmodel35cuts_gender.RData")
save(objects,file="impmodel50cuts_gender.RData")
save(objects,file="impmodel75cuts_gender.RData")

imputed10 = model10$y_imp[[1]]
imputed25 = model25$y_imp[[1]]
imputed35 = model35$y_imp[[1]]
imputed50 = model50$y_imp[[1]]
imputed75 = model75$y_imp[[1]]

write.table(imputed10, "imputed_fbsmodel_gender10.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed25, "imputed_fbsmodel_gender25.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed35, "imputed_fbsmodel_gender35.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed50, "imputed_fbsmodel_gender50.csv", sep = ";", quote = F, col.names=F, row.names=F)
write.table(imputed75, "imputed_fbsmodel_gender75.csv", sep = ";", quote = F, col.names=F, row.names=F)

write.table(model10$bdraw[[1]], "bdraws_fbsmodel10.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model25$bdraw[[1]], "bdraws_fbsmodel25.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model35$bdraw[[1]], "bdraws_fbsmodel35.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model50$bdraw[[1]], "bdraws_fbsmodel50.csv", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model75$bdraw[[1]], "bdraws_fbsmodel75.csv", sep = ";", quote = F, col.names = F, row.names = F)

save(objects,file="impmodel10fbs_gender.RData")
save(objects,file="impmodel25fbs_gender.RData")
save(objects,file="impmodel35fbs_gender.RData")
save(objects,file="impmodel50fbs_gender.RData")
save(objects,file="impmodel75fbs_gender.RData")


----------------------------------------------------------------
library(zoom) 
zm()

setwd("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/bdraws_gender/")
setwd("D:/Documents/downloads/ToreOpsahl/bayes_oxford_random_gender_21052018/bdraws")

bdrawsNULL = read.table("bdraws_cutmodel_nomis.csv", sep = ";", h = F)
bdraws10cut = read.table("bdraws_cutmodel10.csv", sep = ";", h = F)
bdraws25cut = read.table("bdraws_cutmodel25.csv", sep = ";", h = F)
bdraws35cut = read.table("bdraws_cutmodel35.csv", sep = ";", h = F)
bdraws50cut = read.table("bdraws_cutmodel50.csv", sep = ";", h = F)
bdraws75cut = read.table("bdraws_cutmodel75.csv", sep = ";", h = F)
bdraws10fbs = read.table("bdraws_fbsmodel10.csv", sep = ";", h = F)
bdraws25fbs = read.table("bdraws_fbsmodel25.csv", sep = ";", h = F)
bdraws35fbs = read.table("bdraws_fbsmodel35.csv", sep = ";", h = F)
bdraws50fbs = read.table("bdraws_fbsmodel50.csv", sep = ";", h = F)
bdraws75fbs = read.table("bdraws_fbsmodel75.csv", sep = ";", h = F)

plot(density(bdrawsNULL$V1), main = "Intercept")
lines(density(bdraws10cut$V1), col = "red1")
lines(density(bdraws25cut$V1), col = "red2")
lines(density(bdraws35cut$V1), col = "red3")
lines(density(bdraws50cut$V1), col = "red4")
lines(density(bdraws75cut$V1), col = "red4")
lines(density(bdraws10fbs$V1), col = "royalblue")
lines(density(bdraws25fbs$V1), col = "royalblue1")
lines(density(bdraws35fbs$V1), col = "royalblue2")
lines(density(bdraws50fbs$V1), col = "royalblue3")
lines(density(bdraws75fbs$V1), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V2), main = "Popularity")
lines(density(bdraws10cut$V2), col = "red1")
lines(density(bdraws25cut$V2), col = "red2")
lines(density(bdraws35cut$V2), col = "red3")
lines(density(bdraws50cut$V2), col = "red4")
lines(density(bdraws75cut$V2), col = "red4")
lines(density(bdraws10fbs$V2), col = "royalblue")
lines(density(bdraws25fbs$V2), col = "royalblue1")
lines(density(bdraws35fbs$V2), col = "royalblue2")
lines(density(bdraws50fbs$V2), col = "royalblue3")
lines(density(bdraws75fbs$V2), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V3), main = "Year of Study")
lines(density(bdraws10cut$V3), col = "red1")
lines(density(bdraws25cut$V3), col = "red2")
lines(density(bdraws35cut$V3), col = "red3")
lines(density(bdraws50cut$V3), col = "red4")
lines(density(bdraws75cut$V3), col = "red4")
lines(density(bdraws10fbs$V3), col = "royalblue")
lines(density(bdraws25fbs$V3), col = "royalblue1")
lines(density(bdraws35fbs$V3), col = "royalblue2")
lines(density(bdraws50fbs$V3), col = "royalblue3")
lines(density(bdraws75fbs$V3), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V4), main = "Day Active")
lines(density(bdraws10cut$V4), col = "red1")
lines(density(bdraws25cut$V4), col = "red2")
lines(density(bdraws35cut$V4), col = "red3")
lines(density(bdraws50cut$V4), col = "red4")
lines(density(bdraws75cut$V4), col = "red4")
lines(density(bdraws10fbs$V4), col = "royalblue")
lines(density(bdraws25fbs$V4), col = "royalblue1")
lines(density(bdraws35fbs$V4), col = "royalblue2")
lines(density(bdraws50fbs$V4), col = "royalblue3")
lines(density(bdraws75fbs$V4), col = "royalblue3")
zm()

setwd("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/")

# read in the rho estimates for the full data data
effective0 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/rhoestimate_nomissing_cut.txt")

# read in the rho estimates with missing data
effective10cut = read.table("rhoestimates/rhoestimate10cut.txt")
effective10fbs = read.table("rhoestimates/rhoestimate10fbs.txt")

effective25cut = read.table("rhoestimates/rhoestimate25cut.txt")
effective25fbs = read.table("rhoestimates/rhoestimate25fbs.txt")

effective35cut = read.table("rhoestimates/rhoestimate35cut.txt")
effective35fbs = read.table("rhoestimates/rhoestimate35fbs.txt")

effective50cut = read.table("rhoestimates/rhoestimate50cut.txt")
effective50fbs = read.table("rhoestimates/rhoestimate50fbs.txt")

effective75cut = read.table("rhoestimates/rhoestimate75cut.txt")
effective75fbs = read.table("rhoestimates/rhoestimate75fbs.txt")

plot(effective0$V1, type = "l", main = "rho of complete model")
plot(effective10cut$V1, type = "l", main = "rho of CM 10% missing")
plot(effective25cut$V1, type = "l", main = "rho of CM 25% missing")
plot(effective35cut$V1, type = "l", main = "rho of CM 35% missing")
plot(effective50cut$V1, type = "l", main = "rho of CM 50% missing")
plot(effective75cut$V1, type = "l", main = "rho of CM 75% missing")

plot(effective10fbs$V1, type = "l", main = "rho of FBM 10% missing")
plot(effective25fbs$V1, type = "l", main = "rho of FBM 25% missing")
plot(effective35fbs$V1, type = "l", main = "rho of FBM 35% missing")
plot(effective50fbs$V1, type = "l", main = "rho of FBM 50% missing")
plot(effective75fbs$V1, type = "l", main = "rho of FBM 75% missing")

####################################################################################
####################################################################################

imputed10 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_cutmodel_gender10.csv", sep = ";")


z = as.matrix(colMeans(imputed10))
rho = mean(effective10cut$V1[c(1001:16000)])
W10 = W[unique(missing.sample.10$sampleN), unique(missing.sample.10$sampleN)]
X10 = as.matrix(X[unique(missing.sample.10$sampleN),])
b10 = as.matrix(colMeans(bdraws10cut))

# X10 %*% b10
output = as.data.frame(matrix(NA, 1000, length(unique(missing.sample.10$sampleN))))
for(i in 1:1000){
  z = rho * (W10 %*% z)  + X10 %*% b10
  output[i,] = t(z)
  print(i)
}

# maakt dus niets uit, z verandert niet meer.

library(zoom) 
zm()

####################################################################################
# opnieuw runnen en inladen data 

# cut 
imputed10 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_cutmodel_gender10.csv", sep = ";")
imputed25 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_cutmodel_gender25.csv", sep = ";")
imputed35 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_cutmodel_gender35.csv", sep = ";")
imputed50 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_cutmodel_gender50.csv", sep = ";")
imputed75 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_cutmodel_gender75.csv", sep = ";")
# first process this output

# fbs
imputed10 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_fbsmodel_gender10.csv", sep = ";")
imputed25 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_fbsmodel_gender25.csv", sep = ";")
imputed35 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_fbsmodel_gender35.csv", sep = ";")
imputed50 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_fbsmodel_gender50.csv", sep = ";")
imputed75 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/imputed_gender/imputed_fbsmodel_gender75.csv", sep = ";")

imputed = imputed10
imputed = imputed25
imputed = imputed35
imputed = imputed50
imputed = imputed75

# two ways of doing this, take the mean \hat{z} or take the mode indicator function of \hat{z}
mean.z = as.data.frame(as.matrix(colMeans(imputed)))
mean.z$sign <- NA
mean.z[mean.z$V1 < 0,]$sign <- 0
mean.z[mean.z$V1 > 0,]$sign <- 1
table(mean.z$sign)

# geeft dezelfde output als de 4 regels hierboven
#output = matrix(NA, ncol(imputed), 1)
#for(i in 1:ncol(imputed)){
#  person = imputed[,i]
#  person[person < 0] <- 0
#  person[person > 0] <- 1
#  # check
#  
#  if(as.numeric(table(person)[1]) > as.numeric(table(person)[2])){
#    # make it 0
#    output[i,] <- 0
#  } else {
#    # make it 1
#    output[i,] <- 1
#  }
#}

#  combineren van Y en imputed Y
combined.y = as.data.frame(cbind(y[unique(missing.sample.10$sampleN)], mean.z$sign))
combined.y = as.data.frame(cbind(y[unique(missing.sample.25$sampleN)], mean.z$sign))
combined.y = as.data.frame(cbind(y[unique(missing.sample.35$sampleN)], mean.z$sign))
combined.y = as.data.frame(cbind(y[unique(missing.sample.50$sampleN)], mean.z$sign))
combined.y = as.data.frame(cbind(y[unique(missing.sample.75$sampleN)], mean.z$sign))

combined.y$diff = combined.y$V1 - combined.y$V2
sum(abs(combined.y$diff))

#  cut
77  / length(unique(missing.sample.10$sampleN))
191 / length(unique(missing.sample.25$sampleN))
276 / length(unique(missing.sample.35$sampleN))
396 / length(unique(missing.sample.50$sampleN))
686 / length(unique(missing.sample.75$sampleN))

#  cut
77  / length(unique(missing.sample.10$sampleN))
191 / length(unique(missing.sample.25$sampleN))
277 / length(unique(missing.sample.35$sampleN))
399 / length(unique(missing.sample.50$sampleN))
578 / length(unique(missing.sample.75$sampleN))
