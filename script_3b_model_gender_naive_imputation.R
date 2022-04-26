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

# simulate missing data
set.seed(19985)
missing.sample.10 = sample(1899, 185)
set.seed(19985)
missing.sample.25 = sample(1899, 476)
set.seed(19985)
missing.sample.50 = sample(1899, 952)
set.seed(19985)
missing.sample.75 = sample(1899, 1424)

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
"
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
"

########################
### Indirect Mathing ###
########################
yf10 = y; yf10[unique(sort(missing.sample.10))] <- NA
yf25 = y; yf25[unique(sort(missing.sample.25))] <- NA
yf50 = y; yf50[unique(sort(missing.sample.50))] <- NA
yf75 = y; yf75[unique(sort(missing.sample.75))] <- NA

distance.matrix = as.matrix(dist(FBattributes_sel[,-c(8)]))
length(which(distance.matrix[lower.tri(distance.matrix)] == 0)) # are there any persons with equal covariates?

Observed = seq(1,1899)
Observed = Observed[-unique(sort(missing.sample.10))]

#Nobs = dim(distance.matrix)[1]
Nobs = 1899-length(unique(sort(missing.sample.10)))
steps = c(1,3,5,8,13,21,33,54,87)
misclas.table = matrix(NA, Nobs, length(steps))
for(i in 1:Nobs){
  person = Observed[i]
  person.distance = cbind(seq(1,dim(distance.matrix)[1]), distance.matrix[person, ]) # select all distances with person 12 (except person 12)
  person.distance = person.distance[-person,] # remove person 
  person.distance = cbind(person.distance, rank(person.distance[,2])) # rank distances
  person.distance = person.distance[order(person.distance[,3]),]  
  
  misclas_rate = 0
  for(j in 1:length(steps)){
    donors = matrix(person.distance[c(1:steps[j]),], steps[j], 3)
    gender = y[person]
    gender.selected.friends = FBattributes[donors[,1],]$gender01
    overall.friends.gender = ifelse(sum(gender.selected.friends)/length(gender.selected.friends) > .5, 1, 0)
    misclas_rate = ifelse(gender == overall.friends.gender, 0, 1)
    misclas.table[i, j] = misclas_rate
  }
  print(i)
}

misclass = colSums(misclas.table)

output = rbind(cbind(misclas.table[,1], 1),
                cbind(misclas.table[,2], 2),
                cbind(misclas.table[,3], 3),
                cbind(misclas.table[,4], 4),
                cbind(misclas.table[,5], 5),
                cbind(misclas.table[,6], 6),
                cbind(misclas.table[,7], 7),
                cbind(misclas.table[,8], 8),
                cbind(misclas.table[,9], 9))
          
cbind(steps, misclass, misclass/1899)

# So now we take this threshold (only one most similar neighbour), for the 185 cases

Unobserved = sort(unique(missing.sample.10))
output = matrix(NA, 185, 1)
for(i in 1:length(missing.sample.10)){
  person = Unobserved[i]
  person.distance = cbind(seq(1,dim(distance.matrix)[1]), distance.matrix[person, ]) # select all distances with person 12 (except person 12)
  person.distance = person.distance[-person,] # remove person 
  person.distance = cbind(person.distance, rank(person.distance[,2])) # rank distances
  person.distance = person.distance[order(person.distance[,3]),]  
  
  gender = y[person]
  gender.selected.friends = FBattributes[ person.distance[c(1:21), 1],]$gender01
  overall.friends.gender = ifelse(mean(gender.selected.friends) > .5, 1, 0)
  output[i, 1] = ifelse(gender == overall.friends.gender, 0, 1)
}


### old method 
# calculate the best optimal distance to get the best brier score | friends
distance.matrix.missing = distance.matrix[,unique(sort(missing.sample.10))]
distance.matrix.missing = distance.matrix[,unique(sort(missing.sample.25))]
distance.matrix.missing = distance.matrix[,unique(sort(missing.sample.50))]
distance.matrix.missing = distance.matrix[,unique(sort(missing.sample.75))]

missing.sample = missing.sample.10

Nmiss =  length(unique(sort(missing.sample)))
brier.output=matrix(NA,1899-1,Nmiss)
for(i in 1:Nmiss){
  person = unique(sort(missing.sample))[i]
  #friends = unique(LSMI(net, n.seeds = 1, n.neigh = 1, seeds = person, classic = F)$sample)
  person.distance = cbind(seq(1,1898), distance.matrix.missing[-unique(sort(missing.sample))[i],i]) # select all distances with person 12 (except person 12)
  person.distance = cbind(person.distance, rank(person.distance[,2])) # rank distances
  person.distance = person.distance[order(person.distance[,3]),]  
  for(j in 2:1898){
    donors = as.data.frame(person.distance[c(1:j),])
    gender = y[unique(sort(missing.sample))[i]]
    gender.selected.friends = FBattributes[donors[,1],]$gender01
    brier.output[j-1,i] = brier.score(N = 1, f = mean(gender.selected.friends), o = gender)
  }
  print(i)
}

# test print best output
brier.output = brier.output[-1898,]
output = matrix(NA, Nmiss, 1)
for(i in 1:Nmiss){
  output.local = which(brier.output[,i] == min(brier.output[,i]))
  if(length(output.local) == 1){
    output[i,1] = output.local
  } else {
    output[i,1] = output.local[1]
  }
}

hist(output, breaks = 20)

#=> 10 closest persons = 120/185
#=> 20 closest persons = 138/185 

# optimal number of persons >< distance based on brier score < .25
output = matrix(NA, Nmiss, 2)
for(i in 1:Nmiss){
  # output.local = which(brier.output[,i] == min(brier.output[,i])) # take the lowest BS
  #if(length(output.local) == 1){output[i,1] = output.local} else {output[i,1] = output.local[1]}
  
  # take the lowest value untill encountering increasing value
  value = 0
  for(j in 1:1896){
    value=brier.output[j,i]
    if(brier.output[j+1,i] < value){
      value=brier.output[j+1,i]
    } else {
      output[i,1]=j
      output[i,2]=value
      break
    }
  }
}

head(output)
hist(output[,1], main = "Optimal number of donors", xlab = "donors")
hist(output[,2], main = "Brier Score from optimal number of donors", xlab = "Brier Score")
# only 10% of the sample has Brier Score > .5 indicating a small number of persons just cannot 
# be predicted with high accuracy based on the distance
table(ifelse(output[,2]<0.5,0,1))
#write.table(output, "distance.based.brier.score.10.txt", sep = "\t", quote = F, col.names = F, row.names = F)
#write.table(output, "distance.based.brier.score.25.txt", sep = "\t", quote = F, col.names = F, row.names = F)
#write.table(output, "distance.based.brier.score.50.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(output, "distance.based.brier.score.75.txt", sep = "\t", quote = F, col.names = F, row.names = F)

#########################
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

# example
plot(brier.output[,182], main = "example optimal BS = global minimum", xlab = "j", ylab = "Brier Score (BS)")


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
  friends = lsmi(net, n.seed = 1, n.wave = 1, seeds = person)
  friends.data = FBattributes[unlist(friends[[1]][[2]]),]
  #friends.data = friends.data[-person,] # exclude the person we are using to sample
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

########## friends matching optimized ###########
set.seed(19985)
missing.ids = unique(sort(missing.sample.10))
missing.ids = unique(sort(missing.sample.25))
missing.ids = unique(sort(missing.sample.50))
missing.ids = unique(sort(missing.sample.75))

brier.output=matrix(0,length(missing.ids),8)
for(i in 1:length(missing.ids)){
  person = missing.ids[i]
  gender = FBattributes[person,]$gender01
  # lets get the friends of person 1
  for(n in 1:4){
    friends = lsmi(net, n.seed = 1, n.wave = n, seeds = person)
    friends = unique(unlist(friends[[1]][[2]]))
    #friends = friends[-which(friends == person)] # remove person we are studying
    friends.data = FBattributes[friends,] # get attribute data
    friends.gender = mean(friends.data$gender01) # get gender
    brier.output[i,n] = brier.score(N = 1, f = friends.gender, o = gender)
    brier.output[i,n+4] = nrow(friends.data)
    
  }
  print(i)
}

plot(density(brier.output[,5]), main = "amount of friends", xlim = c(0, 1900))
lines(density(brier.output[,6]), col = "blue")
lines(density(brier.output[,7]), col = "red")
lines(density(brier.output[,8]), col = "green")

# how many neighbours is optimal?
output = matrix(NA, length(missing.ids), 2)
for(i in 1:length(missing.ids)){
  optimal = which(brier.output[i,] == min(brier.output[i,]))
  if(length(optimal)>1){optimal = optimal[1]}
  output[i,1] = optimal
  output[i,2] = brier.output[i,optimal]
}

table(output[,1])
write.table(output, "network.based.brier.score.mis10.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(output, "network.based.brier.score.mis25.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(output, "network.based.brier.score.mis50.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(output, "network.based.brier.score.mis75.txt", sep = "\t", quote = F, col.names = F, row.names = F)

################

indirect.brier10 = read.table("distance.based.brier.score.10.txt")
indirect.brier25 = read.table("distance.based.brier.score.25.txt")
indirect.brier50 = read.table("distance.based.brier.score.50.txt")
indirect.brier75 = read.table("distance.based.brier.score.75.txt")

plot(density(indirect.brier10$V2), main = "Brier Score distribution", ylim = c(0, 5), lwd = 3, xlab = "Brier Score")
lines(density(indirect.brier25$V2), col = "blue", lwd = 3)
lines(density(indirect.brier50$V2), col = "red", lwd = 3)
lines(density(indirect.brier75$V2), col = "green", lwd = 3)

##

friend.brier10 = read.table("network.based.brier.score.mis10.txt")
friend.brier25 = read.table("network.based.brier.score.mis25.txt")
friend.brier50 = read.table("network.based.brier.score.mis50.txt")
friend.brier75 = read.table("network.based.brier.score.mis75.txt")

plot(density(friend.brier10$V2), main = "Brier Score distribution", ylim = c(0, 6), xlim = c(0,1), lwd = 3, xlab = "Brier Score")
lines(density(friend.brier25$V2), col = "blue", lwd = 3)
lines(density(friend.brier50$V2), col = "red", lwd = 3)
lines(density(friend.brier75$V2), col = "green", lwd = 3)

table(friend.brier10$V1)
table(friend.brier25$V1)
table(friend.brier50$V1)
table(friend.brier75$V1)




################

