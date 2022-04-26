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

###################
brier.score = function(N=1, f=0.5, o=1){
  # f = probability forecast; o = observed outcome; N = number of forecasts
  brier = 1/N * (f - o)^2
  return(brier)
}

setwd("D:/Documents/downloads/ToreOpsahl")

OClinks_w = read.table("OClinks_w.txt")
colnames(OClinks_w) = c("id1", "id2", "messages")
OClinks_w_g = graph_from_edgelist(as.matrix(OClinks_w[,c(1:2)]), directed = F)
net = igraph_to_network(OClinks_w_g)

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

# MATCH 
# Lets suppose we do not want do perform difficult imputation analyses, 
# and just look at the data to find identical persons in the dataset.

FBattributes$ID = seq(1, nrow(FBattributes))
nrow(FBattributes[!duplicated(FBattributes),]) != nrow(FBattributes)


######################
### Direct Mathing ###
######################
# selecteer identieke cases
library(dplyr)
FBattributes_sel = FBattributes[,-c(1,2,10:15)]
identicals = FBattributes_sel %>% group_by(yearofstudy, 
                                           outdegree,
                                           indegree, 
                                           popular,
                                           averageChar, 
                                           dayactive, 
                                           send.day34th.friends) %>% mutate(n = n()) %>% filter(n > 1) 
identicals = identicals[order(identicals$n),]

identicals = identicals[order(identicals$yearofstudy, 
                              identicals$outdegree,
                              identicals$indegree, 
                              identicals$popular,
                              identicals$averageChar, 
                              identicals$dayactive, 
                              identicals$send.day34th.friends),]


# group the identical cases
group.id = 1
identicals$groupid = 0
for(i in 1:nrow(identicals)){
  if(identicals[i,]$groupid == 0){
    group.amount = identicals[i,]$n
    identicals[i:(i+group.amount-1),]$groupid = group.id
    group.id = group.id + 1
  }
}

# find your person, lets say with missing gender ID = 1360, he/she? has 1 identical person, namely person 1170:
#identicals[identicals$groupid == identicals[identicals$ID == 1360,]$groupid, ]
# check the gender of person 1170:
#FBattributes[FBattributes$ID == 1170 | FBattributes$ID == 1360,]
# one is a female and one is a male, so:
# if we would have taken the gender of person 1170 to be the gender of person 1360 we would have made a mistake

# Now lets check over the entire sample how often we are wrong if we match on identical covariates alone:
Ngroups = max(unique(identicals$groupid))
wrong = 0
for(i in 1:Ngroups){
  group.ids = identicals[identicals$groupid == i, ]$ID
  group.data = FBattributes[group.ids, ]
  if(length(unique(group.data$gender01)) == 2){wrong = wrong+1}
  gender.percentage = sum(group.data$gender01) / nrow(group.data)
}
# in > 51% of the groups we would have made the wrong decision just by matching the identical cases
wrong / Ngroups

########################
### Indirect Mathing ###
########################
brier.output = matrix(0, 1899, 2)
distance.matix = as.matrix(dist(FBattributes_sel[,-c(8)]))
wrong = 0
for(i in 1:1899){
  friends = distance.matix[,i]
  gender = FBattributes[i,]$gender01
  friends = as.data.frame(friends[-i]) # remove person
  friends$ranking = rank(friends)
  # head(friends[order(friends$ranking, decreasing = F),])
  selected.friends = as.numeric(paste(row.names(friends[friends$ranking <= 5,])))
  gender.selected.friends = FBattributes[selected.friends,]$gender01
  if(length(unique(gender.selected.friends)) == 2){wrong = wrong+1}
  brier.output[i,1] = gender - ifelse(mean(gender.selected.friends) < .5, 0, 1)
  brier.output[i,2] = brier.score(N = 1, f = sum(gender.selected.friends) / length(gender.selected.friends), o = FBattributes[i,]$gender01)
}

table(brier.output[,1])
745/1899 # misclassification

# brier score plot
plot(density(brier.output[,2]), 
     main = "density brier score indirect matching", 
     xlab = "Brier Score",
     col = "blue",
     lwd = 5) 

######################
### Friend-Mathing ###
######################
# how about matching based on network structure?
# select a random person from the network to get a missing value
set.seed(19985)
x = c()
brier.output=matrix(0,1899,2)
for(i in 1:1899){
  person = i
  gender = FBattributes[person,]$gender01
  # lets get the friends of person 1
  friends = LSMI(net, n.seeds = 1, n.neigh = 1, seeds = person, classic = F)
  friends.data = FBattributes[unique(friends$sampleN),]
  friends.data = friends.data[-person,] # exclude the person we are using to sample
  # calculate mean gender
  friends.gender = mean(friends.data$gender01)
  brier.output[i,1] = gender - ifelse(friends.gender < .5, 0, 1)
  brier.output[i,2] = brier.score(N = 1, f = friends.gender, o = gender)
}

# misclassification
table(brier.output[,1]) 
1260/1899 # misclassification

# brier score plot
plot(density(brier.output[,2]), 
     main = "density brier score network matching", 
     xlab = "Brier Score",
     col = "blue",
     lwd = 5) 




FBattributes_sel[c(person.distances[c(1:5),2]),]
Ndonors = 10
#FBattributes_sel[c(person.distances[c(1:Ndonors),2]),]$ID
friends.matched = LSMI(net, n.seeds = Ndonors, n.neigh = 1, seeds = FBattributes_sel[c(person.distances[c(1:Ndonors),2]),]$ID, classic = F)

#####
# first find out which persons are more or less similar to the candidate
distance.matrix = as.matrix(dist(FBattributes[,c(2,3,4,5,6,7,8)], diag = TRUE, upper = TRUE))
distance.persons = as.data.frame(t(distance.matrix[candidates,]))

output = list()
for(CAND in 1:length(candidates)){
  ranked.neighbours = row.names(distance.persons[order(distance.persons[,CAND]),])
  ranked.neighbours = ranked.neighbours[ranked.neighbours != candidates[CAND]] # remove the candidate from the list
  
  # find out which of these persons have neighbours that are fairly similar
  # select neighbours of the candidate
  candidate.neighbours = rbind(OClinks_w[OClinks_w$id1==candidates[CAND],],
                               OClinks_w[OClinks_w$id2==candidates[CAND],])           # get the neighbours
  candidate.neighbours = unique(c(candidate.neighbours$id1, candidate.neighbours$id2))# get the neighbours; clean set
  candidate.neighbours = candidate.neighbours[candidate.neighbours!=candidates[CAND]] # remove candidate 

  output[[CAND]] = matrix(NA, distance.cutoff, 5)

  # now select the neighbours of the potential person close to the candidate
  for(POT in 1:distance.cutoff){
    potential.person = ranked.neighbours[POT]
    potential.person.neighbours = rbind(OClinks_w[OClinks_w$id1==potential.person,],
                                        OClinks_w[OClinks_w$id2==potential.person,])                         # get the neighbours
    potential.person.neighbours = unique(c(potential.person.neighbours$id1, potential.person.neighbours$id2))# get the neighbours; clean set
    potential.person.neighbours = potential.person.neighbours[potential.person.neighbours!=potential.person] # remove candidate 
   
    # select the neighbours from the candidate and the neigbours from the potential person
    network.selection = unique(c(candidate.neighbours, potential.person.neighbours))
    network.selection.distance = as.matrix(dist(FBattributes[network.selection ,c(2,3,4,5,6,7,8)]))

    # quantify overlap in neighbours
    overlap = length(network.selection) - length(potential.person.neighbours) + length(candidate.neighbours)
    output[[CAND]][POT,1] = candidates[CAND]
    output[[CAND]][POT,2] = potential.person
    output[[CAND]][POT,3] = overlap
    output[[CAND]][POT,4] = mean(network.selection.distance)
    output[[CAND]][POT,5] = paste(potential.person.neighbours, collapse = " ")
    
  }
  print(CAND)
}

# select the 10 most matching subnetworks
#matching = as.data.frame(output[[1]])
matching = as.data.frame(output[[2]])
matching$V3 = as.numeric(paste(matching$V3))
matching$V4 = as.numeric(paste(matching$V4))
matching = matching[order(matching$V4, decreasing = F) & order(matching$V3, decreasing = T),]
matching = matching[c(1:top),]

# get neighbours in the top 10 matching subnetworks
neighbours = c(as.numeric(paste(as.character(matching$V2))),
               as.numeric(unlist(strsplit(as.character(matching$V5), " "))))
selection = sort(unique(c(candidates[1], neighbours)))
selection = sort(unique(c(candidates[2], neighbours)))

W = W[selection, selection]
FBattributes = FBattributes[selection,]

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

yf_match = y
#yf_match[which(FBattributes$ID == 64)] <- NA
yf_match[which(FBattributes$ID == 519)] <- NA

# X includes YOS 4 without variation (all zero), so remove
X = X[,-c(6)]

#### cut model
source("D:/Documents/networks/scripts/functions/sar_imp_cut.R")
impmodel.cut = sar_combined_mcmc_imp_cut(y = yf_match, x=X, W, ndraw=25000, burn.in=1000, thinning=1,start = start,
                                    prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                    m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                    model = "match")

source("D:/Documents/networks/scripts/functions/sar_imp_fbs.R")
impmodel.fbs = sar_combined_mcmc_imp_fbs(y = yf_match, x=X, W, ndraw = 25000, burn.in = 1000, method = c("probit"),
                                    prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                    m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                    model = "match")

imputed1 = impmodel.cut$y_imp[[1]]
imputed2 = impmodel.fbs$y_imp[[1]]

imputed1[imputed1 <= 0,] <- 0
imputed1[imputed1 > 0,] <- 1
table(imputed1)

imputed2[imputed2 <= 0,] <- 0
imputed2[imputed2 > 0,] <- 1
table(imputed2)





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
                                   model = 10)

model25 = sar_combined_mcmc_imp_fbs(y = yf25, x = X, W = W, ndraw = 25000, burn.in = 1000, method = c("probit"),
                                   prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                   m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                   model = 25)

model35 = sar_combined_mcmc_imp_fbs(y = yf35, x = X, W = W, ndraw = 25000, burn.in = 1000, method = c("probit"),
                                    prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                    m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                    model = 35)

model50 = sar_combined_mcmc_imp_fbs(y = yf50, x = X, W = W, ndraw = 25000, burn.in = 1000, method = c("probit"),
                                    prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0),
                                    m=10, thinning=1, computeMarginalEffects=TRUE, showProgress=TRUE, start = start,
                                    model = 50)

model75 = sar_combined_mcmc_imp_fbs(y = yf75, x = X, W = W, ndraw = 25000, burn.in = 1000, method = c("probit"),
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

# zijn het dezelfde mensen
imputeda = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_cutmodel_gender10.csv", sep = ";")
imputedb = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_fbsmodel_gender10.csv", sep = ";")
# zijn het dezelfde mensen
imputeda = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_cutmodel_gender25.csv", sep = ";")
imputedb = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_fbsmodel_gender25.csv", sep = ";")
# zijn het dezelfde mensen
imputeda = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_cutmodel_gender35.csv", sep = ";")
imputedb = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_fbsmodel_gender35.csv", sep = ";")
# zijn het dezelfde mensen
imputeda = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_cutmodel_gender50.csv", sep = ";")
imputedb = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_fbsmodel_gender50.csv", sep = ";")
# zijn het dezelfde mensen
imputeda = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_cutmodel_gender75.csv", sep = ";")
imputedb = read.table("D:/Documents/downloads/ToreOpsahl/bayes_random_gender_28052018/imputed_random/imputed_fbsmodel_gender75.csv", sep = ";")

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