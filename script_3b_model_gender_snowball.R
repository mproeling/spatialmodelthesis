rm(list = ls(all = TRUE))

# this script is updated by taking Year of Study as a set of dummy variables
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
#set.seed(3030)
#missing.sample.10 = LSMI(net, n.seeds = 12, n.neigh = 1, seeds = NULL, classic = F)
#length(unique(missing.sample.10$sampleN))
#set.seed(6060)
#missing.sample.10 = LSMI(net, n.seeds = 20, n.neigh = 1, seeds = NULL, classic = F)
#length(unique(missing.sample.10$sampleN))

#set.seed(19985)
#missing.sample.25 = LSMI(net, n.seeds = 58, n.neigh = 1, seeds = NULL, classic = F)
#length(unique(missing.sample.25$sampleN))

#set.seed(19985)
#missing.sample.50 = LSMI(net, n.seeds = 159, n.neigh = 1, seeds = NULL, classic = F)
#length(unique(missing.sample.50$sampleN))

#set.seed(19985)
#missing.sample.75 = LSMI(net, n.seeds = 457, n.neigh = 1, seeds = NULL, classic = F)
#length(unique(missing.sample.75$sampleN))

# cluster:
set.seed(19985)
#missing.sample.10 = lsmi(net, n.seed = 11, n.wave = 1, seeds = NULL)
missing.sample.10 = lsmi(net, n.seed = 8, n.wave = 1, seeds = NULL) #edge corrected
length(unique(unlist(missing.sample.10)))

set.seed(19985)
#missing.sample.25 = lsmi(net, n.seed = 44, n.wave = 1, seeds = NULL)
missing.sample.25 = lsmi(net, n.seed = 9, n.wave = 1, seeds = NULL) #edge corrected
length(unique(unlist(missing.sample.25)))

set.seed(19985)
#missing.sample.50 = lsmi(net, n.seed = 159, n.wave = 1, seeds = NULL)
missing.sample.50 = lsmi(net, n.seed = 39, n.wave = 1, seeds = NULL) #edge corrected
length(unique(unlist(missing.sample.50)))

set.seed(19985)
#missing.sample.75 = lsmi(net, n.seed = 457, n.wave = 1, seeds = NULL)
missing.sample.75 = lsmi(net, n.seed = 82, n.wave = 1, seeds = NULL) #edge corrected
length(unique(unlist(missing.sample.75)))

# randomly sample 10 persons with missing data for the clustering method
set.seed(19985) 
missing.sample.cluster = sample(net$n, 10, replace = F)

#######################################
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

# make dummy variables for Year of Study (6 years / categories, so 5 dummy variables, all 0 = year 1)
FBattributes$YOS1 = 0
FBattributes$YOS2 = 0
FBattributes$YOS3 = 0
FBattributes$YOS4 = 0
FBattributes$YOS5 = 0

FBattributes[FBattributes$yearofstudy == 2,]$YOS1 <- 1
FBattributes[FBattributes$yearofstudy == 3,]$YOS2 <- 1
FBattributes[FBattributes$yearofstudy == 4,]$YOS3 <- 1
FBattributes[FBattributes$yearofstudy == 5,]$YOS4 <- 1
FBattributes[FBattributes$yearofstudy == 6,]$YOS5 <- 1

###### read in for testing Gender #####
X = as.matrix(FBattributes[,c(1,
                              which(names(FBattributes) == "indegree"),
                              which(names(FBattributes) == "YOS1"),
                              which(names(FBattributes) == "YOS2"),
                              which(names(FBattributes) == "YOS3"),
                              which(names(FBattributes) == "YOS4"),
                              which(names(FBattributes) == "YOS5"),
                              which(names(FBattributes) == "dayactive"))])

y = FBattributes[,which(names(FBattributes) == "gender01")] 
#y[y > 3] <- 3
y = as.integer(y)

#model1all = glm(y ~ as.matrix(FBattributes[,-c(1,2,10)]), family = "binomial")
#model1gender = glm(y ~ X[,-1], family = "binomial")
## prediction of null model
#newy = predict(model1gender, as.data.frame(X[,-1]), type = "response")
#plot(newy)
#newy[newy > .5] <- 1
#newy[newy <= .5] <- 0
#plot(newy)
#combined.y = y - newy

#source("D:/Documents/networks/scripts/functions/sar_probit.R")
# NULMODEL
#model_null_flag0 = sar_probit_mcmc(y = y, X = X, W = W, ndraw=1000, burn.in=100, showProgress=TRUE,
#                               thinning=1, start = start,
#                               prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0)) # eigen versie met z zelf gemaakt


#source("D:/Documents/networks/scripts/functions/sar_probit_imputenull.R")
#file.remove("D:\\Documents\\downloads\\ToreOpsahl\\rhoestimate_nullmodel.txt")
#model1gender = sar_probit_mcmc(y = y, X = as.matrix(X), W = W, ndraw=25000, burn.in=1000, showProgress=TRUE) # eigen versie met z zelf gemaakt

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

#yf10 = y; yf10[unique(sort(missing.sample.10$sampleN))] <- NA
#yf25 = y; yf25[unique(sort(missing.sample.25$sampleN))] <- NA
#yf50 = y; yf50[unique(sort(missing.sample.50$sampleN))] <- NA
#yf75 = y; yf75[unique(sort(missing.sample.75$sampleN))] <- NA

#cluster:
yf10 = y; yf10[sort(unique(unlist(missing.sample.10)))] <- NA
yf25 = y; yf25[sort(unique(unlist(missing.sample.25)))] <- NA
yf50 = y; yf50[sort(unique(unlist(missing.sample.50)))] <- NA
yf75 = y; yf75[sort(unique(unlist(missing.sample.75)))] <- NA 

yf.cluster = y; yf.cluster[missing.sample.cluster] <- NA

#### cut model
source("D:/Documents/networks/scripts/functions/sar_imp_cut.R")
impmodel.10 = sar_combined_mcmc_imp_cut(y = yf10, x=X, W, ndraw=25000, burn.in=1000, thinning=1,start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = "edge_conditioned_10_flag0")

impmodel.25 = sar_combined_mcmc_imp_cut(y = yf25, x=X, W, ndraw=25000, burn.in=1000, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = "edge_conditioned_25_flag0")

impmodel.50 = sar_combined_mcmc_imp_cut(y = yf50, x=X, W, ndraw=25000, burn.in=1000, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = "edge_conditioned_50_flag0")

impmodel.75 = sar_combined_mcmc_imp_cut(y = yf75, x=X, W, ndraw=25000, burn.in=1000, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = "edge_conditioned_75_flag0")

###
source("D:/Documents/networks/scripts/functions/sar_imp_cut_donor_new.R")
impmodel.cluster.donor = sar_combined_mcmc_imp_cut_cluster(y = yf.cluster, x=X, W, ndraw=1000, burn.in=100, thinning=1,start = start,
                                                           prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                                           m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                                           model = "cluster10", Nd = 5)
write.table(impmodel.cluster.donor$y_imp[[1]], "y_imp_cluster.csv", sep = ";", quote = F, col.names = F, row.names = F)



impmodel.cluster.normal.10 = sar_combined_mcmc_imp_cut(y = yf.cluster, x=X, W, ndraw=10000, burn.in=100, thinning=1,start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = "cluster10.normal")

write.table(impmodel.cluster.normal.10$y_imp[[1]], "y_imp_normalcut.csv", sep = ";", quote = F, col.names = F, row.names = F)

impmodel.cluster.normal.fbs.10 = sar_combined_mcmc_imp_fbs(y = yf.cluster, x=X, W, ndraw=10000, burn.in=100, thinning=1,start = start,
                                                       prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                                       m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                                       model = "cluster10.normal.fbs")

write.table(impmodel.cluster.normal.fbs.10$y_imp[[1]], "y_imp_normalfbs.csv", sep = ";", quote = F, col.names = F, row.names = F)




####################
output = read.table("clusteroutput.txt", sep = "\t")
output.new = matrix(0, 5000, 3)
count = 1
for(i in seq(1,nrow(output),3)){
  output.new[count,1] = as.numeric(paste(output[i,]))
  output.new[count,2] = as.numeric(paste(output[i+1,]))
  output.new[count,3] = paste(output[i+2,])
  count = count + 1
  print(i)
}
output.new = as.data.frame(output.new)

first = output.new[output.new$V1 == 1,]
second = output.new[output.new$V1 == 2,]
third = output.new[output.new$V1 == 3,]
fourth = output.new[output.new$V1 == 4,]
fifth = output.new[output.new$V1 == 5,]
sixth = output.new[output.new$V1 == 6,]
seven = output.new[output.new$V1 == 7,]
eight = output.new[output.new$V1 == 8,]
nine = output.new[output.new$V1 == 9,]
ten = output.new[output.new$V1 == 10,]

table(is.na(first$V2))
table(is.na(second$V2))
table(is.na(third$V2))
table(is.na(fourth$V2))
table(is.na(fifth$V2))
table(is.na(sixth$V2))
table(is.na(seven$V2))
table(is.na(eight$V2))
table(is.na(nine$V2))
table(is.na(ten$V2))

####################
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
model10 = sar_combined_mcmc_imp_fbs(y = yf10, x = X, W = W, ndraw = 25000, burn.in = 1000, method = c("probit"),
                                   prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                   m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                   model = "edge_conditioned_10_flag")

model25 = sar_combined_mcmc_imp_fbs(y = yf25, x = X, W = W, ndraw = 25000, burn.in = 1000, method = c("probit"),
                                   prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                   m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                   model = "edge_conditioned_25_flag0")

model50 = sar_combined_mcmc_imp_fbs(y = yf50, x = X, W = W, ndraw = 25000, burn.in = 1000, method = c("probit"),
                                    prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                    m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                    model = "edge_conditioned_50_flag0")

model75 = sar_combined_mcmc_imp_fbs(y = yf75, x = X, W = W, ndraw = 25000, burn.in = 1000, method = c("probit"),
                                    prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                    m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                    model = "edge_conditioned_75_flag0")


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

setwd("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/bdraws_snowball//")

bdrawsNULL = read.table("bdraws_cutmodel_nomis.csv", sep = ";", h = F)
bdraws10cut = read.table("bdraws_cutmodel10.csv", sep = ";", h = F)
bdraws25cut = read.table("bdraws_cutmodel25.csv", sep = ";", h = F)
#bdraws35cut = read.table("bdraws_cutmodel35.csv", sep = ";", h = F)
bdraws50cut = read.table("bdraws_cutmodel50.csv", sep = ";", h = F)
bdraws75cut = read.table("bdraws_cutmodel75.csv", sep = ";", h = F)
bdraws10fbs = read.table("bdraws_fbsmodel10.csv", sep = ";", h = F)
bdraws25fbs = read.table("bdraws_fbsmodel25.csv", sep = ";", h = F)
#bdraws35fbs = read.table("bdraws_fbsmodel35.csv", sep = ";", h = F)
bdraws50fbs = read.table("bdraws_fbsmodel50.csv", sep = ";", h = F)
bdraws75fbs = read.table("bdraws_fbsmodel75.csv", sep = ";", h = F)

plot(density(bdrawsNULL$V1), main = "Intercept", lwd = 3)
lines(density(bdraws10cut$V1), col = "red1")
lines(density(bdraws25cut$V1), col = "red2")
#lines(density(bdraws35cut$V1), col = "red3")
lines(density(bdraws50cut$V1), col = "red3")
lines(density(bdraws75cut$V1), col = "red4")
lines(density(bdraws10fbs$V1), col = "royalblue")
lines(density(bdraws25fbs$V1), col = "royalblue1")
#lines(density(bdraws35fbs$V1), col = "royalblue2")
lines(density(bdraws50fbs$V1), col = "royalblue2")
lines(density(bdraws75fbs$V1), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V2), main = "Indegree",  lwd = 3)
lines(density(bdraws10cut$V2), col = "red1")
lines(density(bdraws25cut$V2), col = "red2")
#lines(density(bdraws35cut$V2), col = "red3")
lines(density(bdraws50cut$V2), col = "red4")
lines(density(bdraws75cut$V2), col = "red4")
lines(density(bdraws10fbs$V2), col = "royalblue")
lines(density(bdraws25fbs$V2), col = "royalblue1")
#lines(density(bdraws35fbs$V2), col = "royalblue2")
lines(density(bdraws50fbs$V2), col = "royalblue3")
lines(density(bdraws75fbs$V2), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V3), main = "Year of Study 2", lwd = 3)
lines(density(bdraws10cut$V3), col = "red1")
lines(density(bdraws25cut$V3), col = "red2")
#lines(density(bdraws35cut$V3), col = "red3")
lines(density(bdraws50cut$V3), col = "red3")
lines(density(bdraws75cut$V3), col = "red4")
lines(density(bdraws10fbs$V3), col = "royalblue")
lines(density(bdraws25fbs$V3), col = "royalblue1")
#lines(density(bdraws35fbs$V3), col = "royalblue2")
lines(density(bdraws50fbs$V3), col = "royalblue2")
lines(density(bdraws75fbs$V3), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V4), main = "Year of Study 3", lwd = 3)
lines(density(bdraws10cut$V4), col = "red1")
lines(density(bdraws25cut$V4), col = "red2")
#lines(density(bdraws35cut$V4), col = "red3")
lines(density(bdraws50cut$V4), col = "red3")
lines(density(bdraws75cut$V4), col = "red4")
lines(density(bdraws10fbs$V4), col = "royalblue")
lines(density(bdraws25fbs$V4), col = "royalblue1")
#lines(density(bdraws35fbs$V4), col = "royalblue2")
lines(density(bdraws50fbs$V4), col = "royalblue2")
lines(density(bdraws75fbs$V4), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V5), main = "Year of Study 4", lwd = 3)
lines(density(bdraws10cut$V5), col = "red1")
lines(density(bdraws25cut$V5), col = "red2")
#lines(density(bdraws35cut$V5), col = "red3")
lines(density(bdraws50cut$V5), col = "red3")
lines(density(bdraws75cut$V5), col = "red4")
lines(density(bdraws10fbs$V5), col = "royalblue")
lines(density(bdraws25fbs$V5), col = "royalblue1")
#lines(density(bdraws35fbs$V5), col = "royalblue2")
lines(density(bdraws50fbs$V5), col = "royalblue2")
lines(density(bdraws75fbs$V5), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V6), main = "Year of Study 5", lwd = 3)
lines(density(bdraws10cut$V6), col = "red1")
lines(density(bdraws25cut$V6), col = "red2")
#lines(density(bdraws35cut$V6), col = "red3")
lines(density(bdraws50cut$V6), col = "red3")
lines(density(bdraws75cut$V6), col = "red4")
lines(density(bdraws10fbs$V6), col = "royalblue")
lines(density(bdraws25fbs$V6), col = "royalblue1")
#lines(density(bdraws35fbs$V6), col = "royalblue2")
lines(density(bdraws50fbs$V6), col = "royalblue2")
lines(density(bdraws75fbs$V6), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V7), main = "Year of Study 6", lwd = 3)
lines(density(bdraws10cut$V7), col = "red1")
lines(density(bdraws25cut$V7), col = "red2")
#lines(density(bdraws35cut$V7), col = "red3")
lines(density(bdraws50cut$V7), col = "red3")
lines(density(bdraws75cut$V7), col = "red4")
lines(density(bdraws10fbs$V7), col = "royalblue")
lines(density(bdraws25fbs$V7), col = "royalblue1")
#lines(density(bdraws35fbs$V7), col = "royalblue2")
lines(density(bdraws50fbs$V7), col = "royalblue2")
lines(density(bdraws75fbs$V7), col = "royalblue3")
zm()

plot(density(bdrawsNULL$V8), main = "Day Active", lwd = 3)
lines(density(bdraws10cut$V8), col = "red1")
lines(density(bdraws25cut$V8), col = "red2")
#lines(density(bdraws35cut$V8), col = "red3")
lines(density(bdraws50cut$V8), col = "red3")
lines(density(bdraws75cut$V8), col = "red4")
lines(density(bdraws10fbs$V8), col = "royalblue")
lines(density(bdraws25fbs$V8), col = "royalblue1")
#lines(density(bdraws35fbs$V8), col = "royalblue2")
lines(density(bdraws50fbs$V8), col = "royalblue2")
lines(density(bdraws75fbs$V8), col = "royalblue3")
zm()

# read in the rho estimates for the full data data
setwd("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/")
                        
effective0 = read.table("rhoestimatesnowball/rhoestimate_nomissing_cut.txt")

# read in the rho estimates with missing data
effective10cut = read.table("rhoestimatesnowball/rhoestimate10cut.txt")
effective10fbs = read.table("rhoestimatesnowball/rhoestimate10fbs.txt")

effective25cut = read.table("rhoestimatesnowball/rhoestimate25cut.txt")
effective25fbs = read.table("rhoestimatesnowball/rhoestimate25fbs.txt")

#effective35cut = read.table("rhoestimatesnowball/rhoestimate35cut.txt")
#effective35fbs = read.table("rhoestimatesnowball/rhoestimate35fbs.txt")

effective50cut = read.table("rhoestimatesnowball/rhoestimate50cut.txt")
effective50fbs = read.table("rhoestimatesnowball/rhoestimate50fbs.txt")

effective75cut = read.table("rhoestimatesnowball/rhoestimate75cut.txt")
effective75fbs = read.table("rhoestimatesnowball/rhoestimate75fbs.txt")

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

imputed10 = read.table("imputed_snowball/imputed_cutmodel_gender10.csv", sep = ";")


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

imputednull = read.table("D:/Documents/downloads/ToreOpsahl/imputednonmissing/imputednonmissing.csv", sep = ";")

# cut 
imputed10 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_cutmodel_gender10.csv", sep = ";")
imputed25 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_cutmodel_gender25.csv", sep = ";")
imputed35 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_cutmodel_gender35.csv", sep = ";")
imputed50 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_cutmodel_gender50.csv", sep = ";")
imputed75 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_cutmodel_gender75.csv", sep = ";")
# first process this output

# fbs
imputed10 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_fbsmodel_gender10.csv", sep = ";")
imputed25 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_fbsmodel_gender25.csv", sep = ";")
imputed35 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_fbsmodel_gender35.csv", sep = ";")
imputed50 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_fbsmodel_gender50.csv", sep = ";")
imputed75 = read.table("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018/imputed_snowball/imputed_fbsmodel_gender75.csv", sep = ";")

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
combined.y = as.data.frame(cbind(y, mean.z$sign))
combined.y = as.data.frame(cbind(y[unique(missing.sample.10$sampleN)], mean.z$sign))
combined.y = as.data.frame(cbind(y[unique(missing.sample.25$sampleN)], mean.z$sign))
combined.y = as.data.frame(cbind(y[unique(missing.sample.35$sampleN)], mean.z$sign))
combined.y = as.data.frame(cbind(y[unique(missing.sample.50$sampleN)], mean.z$sign))
combined.y = as.data.frame(cbind(y[unique(missing.sample.75$sampleN)], mean.z$sign))

combined.y = as.data.frame(cbind(y[sort(unique(unlist(missing.sample.10)))], mean.z$sign)) # cluster

combined.y = as.data.frame(cbind(y[sort(missing.sample.10)], mean.z$sign))


combined.y$diff = combined.y$y - combined.y$V2
combined.y$diff = combined.y$V1 - combined.y$V2
sum(abs(combined.y$diff))

#  nonmissing
sum(abs(combined.y$diff))
629 / length(y)

#  cut
113  / length(unique(missing.sample.10$sampleN))
246 / length(unique(missing.sample.25$sampleN))
371 / length(unique(missing.sample.35$sampleN))
410 / length(unique(missing.sample.50$sampleN))
653 / length(unique(missing.sample.75$sampleN))

#  fbs
115  / length(unique(missing.sample.10$sampleN))
232 / length(unique(missing.sample.25$sampleN))
363 / length(unique(missing.sample.35$sampleN))
421 / length(unique(missing.sample.50$sampleN))
591 / length(unique(missing.sample.75$sampleN))

576 / 1411

# two ways of doing this, take the mean \hat{z} or take the mode indicator function of \hat{z}
mean.z = as.data.frame(as.matrix(colMeans(imputeda)))
mean.z$sign <- NA
mean.z[mean.z$V1 < 0,]$sign <- 0
mean.z[mean.z$V1 > 0,]$sign <- 1
table(mean.z$sign)

combined.y1 = as.data.frame(cbind(y[unique(missing.sample.10$sampleN)], mean.z$sign))
combined.y1 = as.data.frame(cbind(y[unique(missing.sample.25$sampleN)], mean.z$sign))
combined.y1 = as.data.frame(cbind(y[unique(missing.sample.35$sampleN)], mean.z$sign))
combined.y1 = as.data.frame(cbind(y[unique(missing.sample.50$sampleN)], mean.z$sign))
combined.y1 = as.data.frame(cbind(y[unique(missing.sample.75$sampleN)], mean.z$sign))

# two ways of doing this, take the mean \hat{z} or take the mode indicator function of \hat{z}
mean.z = as.data.frame(as.matrix(colMeans(imputedb)))
mean.z$sign <- NA
mean.z[mean.z$V1 < 0,]$sign <- 0
mean.z[mean.z$V1 > 0,]$sign <- 1
table(mean.z$sign)

combined.y2 = as.data.frame(cbind(y[unique(missing.sample.10$sampleN)], mean.z$sign))
combined.y2 = as.data.frame(cbind(y[unique(missing.sample.25$sampleN)], mean.z$sign))
combined.y2 = as.data.frame(cbind(y[unique(missing.sample.35$sampleN)], mean.z$sign))
combined.y2 = as.data.frame(cbind(y[unique(missing.sample.50$sampleN)], mean.z$sign))
combined.y2 = as.data.frame(cbind(y[unique(missing.sample.75$sampleN)], mean.z$sign))

test = cbind(combined.y1, combined.y2)
test = test[,-3]
table(test[,2] - test[,3])

################## analyse clustered output #############
cluster.imp = read.table("y_imp_cluster.csv", sep = ";")
cut.imp = read.table("y_imp_normalcut.csv", sep = ";")

setwd("D:/Documents/downloads/ToreOpsahl/")
clusteroutput = read.table("clusteroutput.txt", sep = "\t")
colnames(clusteroutput) = c("personid", "N1person", "donorid", "N1donors", "Nshared", 
                            "ratiosharedIsharedD", "rationotsharedInotsharedD", "distance", "selected")

persons = unique(clusteroutput$personid)
for(i in 1:10){
  person1 = clusteroutput[clusteroutput$personid == persons[i],]  
  sort(table(person1$donorid))
}
# select 1rst person


##### bias
imputed = impmodel.10$y_imp[[1]]
imputed = impmodel.25$y_imp[[1]]
imputed = impmodel.50$y_imp[[1]]
imputed = impmodel.75$y_imp[[1]]

imputed = model10$y_imp[[1]]
imputed = model25$y_imp[[1]]
imputed = model50$y_imp[[1]]
imputed = model75$y_imp[[1]]

mean.z = as.data.frame(as.matrix(colMeans(imputed)))
mean.z$sign <- NA
mean.z[mean.z$V1 < 0,]$sign <- 0
mean.z[mean.z$V1 > 0,]$sign <- 1
table(mean.z$sign)

# combined.y = as.data.frame(cbind(y[sort(unique(unlist(missing.sample.10)))], mean.z$sign))
combined.y = as.data.frame(cbind(sort(unique(unlist(missing.sample.75))),
                                 y[sort(unique(unlist(missing.sample.75)))], 
                                 mean.z$sign))
combined.y$diff = combined.y$V2 - combined.y$V3
sum(abs(combined.y$diff)) / nrow(combined.y)
write.table(combined.y, "combinedy.cut75.edgecorrected.csv", sep = ";", quote = F, col.names = F, row.names = F)
## save bdraws

write.table(impmodel.10$bdraw[[1]], "bdraw_snow_10_flag2", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.25$bdraw[[1]], "bdraw_snow_25_flag2", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.50$bdraw[[1]], "bdraw_snow_50_flag1", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.75$bdraw[[1]], "bdraw_snow_75_flag2", sep = ";", quote = F, col.names = F, row.names = F)

write.table(model10$bdraw[[1]], "bdraw_snow_10_fbs_flag2", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model25$bdraw[[1]], "bdraw_snow_25_fbs_flag0", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model50$bdraw[[1]], "bdraw_snow_50_fbs_flag0", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model75$bdraw[[1]], "bdraw_snow_75_fbs_flag0", sep = ";", quote = F, col.names = F, row.names = F)


#######
setwd("D:/Documents/downloads/ToreOpsahl/bayes_snowball_flagcomparison/")
cut = read.table("combinedy.cut10.csv", sep = ";")
fbs = read.table("combinedy.fbs10.csv", sep = ";")

cut = read.table("combinedy.cut25.csv", sep = ";")
fbs = read.table("combinedy.fbs25.csv", sep = ";")

cut = read.table("combinedy.cut50.csv", sep = ";")
fbs = read.table("combinedy.fbs50.csv", sep = ";")

cut = read.table("combinedy.cut75.csv", sep = ";")
fbs = read.table("combinedy.fbs75.csv", sep = ";")

test = cbind(cut, fbs)
table(test[,2] - test[,6])  # gender everybody identical?
table(test[,4], test[,8])   # misclassification in off diagonal
test$diff = test[,4]^2 + test[,8]^2
#test[test$diff > 1, ]$diff <- 1 # recode all to 1
table(test$diff)
table(test[,2], test$diff)
# is one model predicting better compared to the other model?
table(test[,2], test[,4]) 
table(test[,6], test[,8]) 


setwd("D:/Documents/downloads/ToreOpsahl/bayes_random_flagcomparison/")

test = read.table("")

write.table(impmodel.10$y_imp[[1]], "yimp_snow_10_cut_flag2", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.25$y_imp[[1]], "yimp_snow_25_cut_flag2", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.50$y_imp[[1]], "yimp_snow_50_cut_flag2", sep = ";", quote = F, col.names = F, row.names = F)
write.table(impmodel.75$y_imp[[1]], "yimp_snow_75_cut_flag2", sep = ";", quote = F, col.names = F, row.names = F)

write.table(model10$y_imp[[1]], "yimp_snow_10_fbs_flag2", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model25$y_imp[[1]], "yimp_snow_25_fbs_flag2", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model50$y_imp[[1]], "yimp_snow_50_fbs_flag2", sep = ";", quote = F, col.names = F, row.names = F)
write.table(model75$y_imp[[1]], "yimp_snow_75_fbs_flag2", sep = ";", quote = F, col.names = F, row.names = F)



###################
brier.score = function(x, N=1, f=0.5, o=1){
  # f = probability forecast; o = observed outcome; N = number of forecasts
  brier = 1/N * (x$f - x$o)^2
  return(brier)
}

#impmodel.10 = read.table("yimp_snow_75_fbs_flag2", sep = ";")
impmodel = impmodel.75$y_imp[[1]]
z.indicator = ifelse(impmodel > 0, 1, 0)
z.perc = colSums(z.indicator) / nrow(z.indicator)
z.table = cbind(z.perc, y[unique(unlist(missing.sample.75))])
colnames(z.table) = c("f", "o")
write.table(brier.score(as.data.frame(z.table)), "brier.score.75snow_edge_cut_flag0.txt", quote = F, sep = ";", row.names = F, col.names=F)
write.table(brier.score(as.data.frame(z.table)), "brier.score.10snow_edge_fbs_flag0.txt", quote = F, sep = ";", row.names = F, col.names=F)
##
rm(list = ls(all = TRUE))
setwd("D:\\Documents\\downloads\\ToreOpsahl\\bayes_random_flagcomparison")

brier.values10 = read.table("brier.score.10random_flag0.txt", sep = ";", h = F)
brier.values25 = read.table("brier.score.25random_flag0.txt", sep = ";", h = F)
brier.values50 = read.table("brier.score.50random_flag0.txt", sep = ";", h = F)
brier.values75 = read.table("brier.score.75random_flag0.txt", sep = ";", h = F)

brier.values10 = read.table("brier.score.fbs.10random_flag0.txt", sep = ";", h = F)
brier.values25 = read.table("brier.score.fbs.25random_flag0.txt", sep = ";", h = F)
brier.values50 = read.table("brier.score.fbs.50random_flag0.txt", sep = ";", h = F)
brier.values75 = read.table("brier.score.fbs.75random_flag0.txt", sep = ";", h = F)

mean(brier.values10$V1)
mean(brier.values25$V1)
mean(brier.values50$V1)
mean(brier.values75$V1)


setwd("D:\\Documents\\downloads\\ToreOpsahl\\bayes_random_flagcomparison")
brier.values10 = read.table("brier.score.10snow_flag0.txt", sep = ";", h = F)
brier.values25 = read.table("brier.score.25snow_flag0.txt", sep = ";", h = F)
brier.values50 = read.table("brier.score.50snow_flag0.txt", sep = ";", h = F)
brier.values75 = read.table("brier.score.75snow_flag0.txt", sep = ";", h = F)

brier.values10 = read.table("brier.score.10snow_fbs_flag0.txt", sep = ";", h = F)
brier.values25 = read.table("brier.score.25snow_fbs_flag0.txt", sep = ";", h = F)
brier.values50 = read.table("brier.score.50snow_fbs_flag0.txt", sep = ";", h = F)
brier.values75 = read.table("brier.score.75snow_fbs_flag0.txt", sep = ";", h = F)

setwd("D:\\Documents\\downloads\\ToreOpsahl\\bayes_snowball_edgecorrected")
brier.values10 = read.table("brier.score.10snow_edge_cut_flag0.txt", sep = ";", h = F)
brier.values25 = read.table("brier.score.25snow_edge_cut_flag0.txt", sep = ";", h = F)
brier.values50 = read.table("brier.score.50snow_edge_cut_flag0.txt", sep = ";", h = F)
brier.values75 = read.table("brier.score.75snow_edge_cut_flag0.txt", sep = ";", h = F)

brier.values10 = read.table("brier.score.10snow_edge_fbs_flag0.txt", sep = ";", h = F)
brier.values25 = read.table("brier.score.25snow_edge_fbs_flag0.txt", sep = ";", h = F)
brier.values50 = read.table("brier.score.50snow_edge_fbs_flag0.txt", sep = ";", h = F)
brier.values75 = read.table("brier.score.75snow_edge_fbs_flag0.txt", sep = ";", h = F)

plot(density(brier.values10$V1), 
     main = "full bayes with grid approximation",
     xlab = "brier score", ylim = c(0, 4), lwd = 5)
lines(density(brier.values25$V1), col = "blue", lwd = 5)
lines(density(brier.values50$V1), col = "green", lwd = 5) 
lines(density(brier.values75$V1), col = "red", lwd = 5) 
               
main = "cut model with grid evaluation",
main = "full bayes with spline approximation",
main = "full bayes with Chebyshev approximation",




# read in the rho estimates with missing data
effective10cut = read.table("rhoestimateedge_conditioned_10_flag0cut.txt")
effective10fbs = read.table("rhoestimateedge_conditioned_10_flag0fbs.txt")

effective25cut = read.table("rhoestimateedge_conditioned_25_flag0cut.txt")
effective25fbs = read.table("rhoestimateedge_conditioned_25_flag0fbs.txt")

effective50cut = read.table("rhoestimateedge_conditioned_50_flag0cut.txt")
effective50fbs = read.table("rhoestimateedge_conditioned_50_flag0fbs.txt")

effective75cut = read.table("rhoestimateedge_conditioned_75_flag0cut.txt")
effective75fbs = read.table("rhoestimateedge_conditioned_75_flag0fbs.txt")

plot(effective10cut$V1, type = "l", main = "rho of CM 10% missing", xlab = "draws", ylab = "rho")
plot(effective25cut$V1, type = "l", main = "rho of CM 25% missing", xlab = "draws", ylab = "rho")
plot(effective50cut$V1, type = "l", main = "rho of CM 50% missing", xlab = "draws", ylab = "rho")
plot(effective75cut$V1, type = "l", main = "rho of CM 75% missing", xlab = "draws", ylab = "rho")

plot(effective10fbs$V1, type = "l", main = "rho of FBM 10% missing", xlab = "draws", ylab = "rho")
plot(effective25fbs$V1, type = "l", main = "rho of FBM 25% missing", xlab = "draws", ylab = "rho")
plot(effective50fbs$V1, type = "l", main = "rho of FBM 50% missing", xlab = "draws", ylab = "rho")
plot(effective75fbs$V1, type = "l", main = "rho of FBM 75% missing", xlab = "draws", ylab = "rho")


