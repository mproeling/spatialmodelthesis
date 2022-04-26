rm(list = ls(all = TRUE))

library(igraph)


setwd("D:/Documents/downloads/ToreOpsahl/")

# data from:
# https://toreopsahl.com/datasets/#online_social_network
# Weighted longitudinal one-mode network (weighted by number of characters): OCnodeslinks_chars.txt
# Binary longitudinal one-mode network: OCnodeslinks.txt
# Weighted static one-mode network (weighted by number of characters): OClinks_w_chars.txt
# Weighted static one-mode network (weighted by number of messages): OClinks_w.txt


OClinks_w = read.table("OClinks_w.txt")
colnames(OClinks_w) = c("id1", "id2", "messages")

W = get.adjacency(graph.edgelist(as.matrix(OClinks_w[,c(1:2)]),directed=TRUE), sparse = igraph_opt("sparsematrices"))
# write assymmetric matrix 
writeMM(messages.W, file='Wmatrix.01.txt')
# make the matrix W symmetric 
symmetric01.W = W
# W must be symmetrix, so we take the highest number (messages received or send, to the other person)
for(i in 1:dim(symmetric01.W)[1]){ # rows
  for(j in 1:dim(symmetric01.W)[2]){ # columns
    if(j > i){
      send = symmetric01.W[i,j]
      received = symmetric01.W[j,i]
      symmetric01.W[i,j] = max(send, received)
      symmetric01.W[j,i] = max(send, received)
    }
  }
  print(i)
}

writeMM(symmetric01.W, file='Wmatrix.01.sym.txt')

# make row stochastic
sumofrows = rowSums(symmetric01.W)
for(i in 1:dim(symmetric01.W)[1]){
  if(sumofrows[i] != 1){
    value = 1 / sumofrows[i]
    for(j in 1:dim(symmetric01.W)[1]){
      if(symmetric01.W[i,j] != 0) symmetric01.W[i,j] = value 
    }
  }
  print(i)
}
writeMM(symmetric01.W, file='Wmatrix.01.sym.rowstochastic.txt')

####################################################
# add weights amount of messages send and received
messages.W = W
for(i in 1:dim(messages.W)[1]){
  for(j in 1:dim(messages.W)[2]){
    if(j > i){
      # get weights for upper (select1) and lower (select2) triangular
      select1 = OClinks_w[OClinks_w$id1 == i,]
      select1 = select1[select1$id2 == j,]
      select2 = OClinks_w[OClinks_w$id1 == j,]
      select2 = select2[select2$id2 == i,]
      
      # combine weights and take the mean
      select = rbind(select1, select2)
      if(nrow(select) != 0){
        select$messages = mean(select$messages)
      }
      
      if(nrow(select1) != 0) messages.W[i,j] = unique(select$messages)
      if(nrow(select2) != 0) messages.W[j,i] = unique(select$messages)
    }
  }
  print(i)
}

writeMM(messages.W, file='Wmatrix.meanmessages.sym.txt')
messages.W[c(1:10), c(1:10)] # test

#messages.W = as(as.matrix(readMM(file='Wmatrix.meanmessages.sym.txt')), "dgCMatrix") 
# make row stochastic
sumofrows = rowSums(messages.W)
stoch.messages.W = messages.W
for(i in 1:dim(messages.W)[1]){
  if(sumofrows[i] > 1){
    # start to loop over the values of i
    for(j in 1:dim(messages.W)[1]){
      value = stoch.messages.W[i,j]
      if(value != 0){
        total = sumofrows[i]
        stoch.value = value / total
        stoch.messages.W[i,j] = stoch.value 
      }
    }
  }
  print(i)
}

stoch.messages.W[c(1:10), c(1:10)] # test
table(rowSums(stoch.messages.W)) # check 
writeMM(stoch.messages.W, file='Wmatrix.meanmessages.rowstochastic.txt')


####


####################################################
# add weights amount of messages send and received MAX
messages.W = W
for(i in 1:dim(messages.W)[1]){
  for(j in 1:dim(messages.W)[2]){
    if(j > i){
      # get weights for upper (select1) and lower (select2) triangular
      select1 = OClinks_w[OClinks_w$id1 == i,]
      select1 = select1[select1$id2 == j,]
      select2 = OClinks_w[OClinks_w$id1 == j,]
      select2 = select2[select2$id2 == i,]
      
      # combine weights and take the mean
      select = rbind(select1, select2)
      if(nrow(select) != 0){
        select$messages = max(select$messages)
      }
      
      if(nrow(select1) != 0) messages.W[i,j] = unique(select$messages)
      if(nrow(select2) != 0) messages.W[j,i] = unique(select$messages)
    }
  }
  print(i)
}

writeMM(messages.W, file='Wmatrix.meanmessages.sym.txt')
messages.W[c(1:10), c(1:10)] # test



# write.table(messages.W, "Wmatrix.messages.txt", sep = ";", quote = F, col.names = F, row.names = F)

attributes = read.table("OCattrib_1899_genderYearofstudy.txt", h = T)
attributes$outdegree = rowSums(messages.W)
attributes$indegree = colSums(messages.W)
attributes$popular = attributes$indegree + attributes$outdegree

# there are 549 people with outdegree 0 and 37 with indegree 0
nrow(attributes[attributes$outdegree == 0,])
nrow(attributes[attributes$indegree == 0,])
nrow(attributes[attributes$outdegree == 0 & attributes$indegree == 0,]) # but there are no persons with both 0
# this means that 549 recieved a message but did not respond, and 37 send a message but did not receive a response. 

# W must be symmetrix, so we take the highest number (messages received or send, to the other person)
for(i in 1:dim(messages.W)[1]){ # rows
  for(j in 1:dim(messages.W)[2]){ # columns
    if(j > i){
      send = messages.W[i,j]
      received = messages.W[j,i]
      messages.W[i,j] = mean(send, received)
      messages.W[j,i] = mean(send, received)
    }
  }
  print(i)
}

writeMM(messages.W, file='Wmatrix.messages.sym.txt')

# we cant use chars and messages in one model because chars will only be >0 for persons 
# who received or sent a message. 
# OClinks_c = read.table("OClinks_w_chars.txt")
OClinks_c = read.table("OCnodeslinks_chars.txt")
colnames(OClinks_c) = c("time", "id1", "id2", "characters")
# Self-loops in the longitudinal edgelist signal the time that users registered on the site.
# set to date (remove times)
OClinks_c$datetime = as.POSIXct(OClinks_c$time, format = "%Y-%m-%d %H:%M:%S")
OClinks_c$time = as.Date(OClinks_c$time)

# order file on time
OClinks_c = OClinks_c[order(OClinks_c$time),]
OClinks_c$datetime[nrow(OClinks_c)] - OClinks_c$datetime[1]

# 3-week lifespan because 3 weeks represent the best approximation of the time at which 
# the rates of increase in messages and in new tie formation stabilize, while the system 
# is still rapidly growing
i = 1
OClinks_c$order = seq(1, nrow(OClinks_c))
OClinks_c$window = 0
for(j in 2:nrow(OClinks_c)){
  day1 = round(OClinks_c$time[i]) 
  day2 = round(OClinks_c$time[j]) 
  window = difftime(day2, day1, units = "days"); window
  print(j)
  if(window > 3*7){
    window.data = OClinks_c[i:j,]
    # add other days with the same date as the end date
    window.data = unique(rbind(window.data, OClinks_c[OClinks_c$time == tail(window.data, n = 1)$time, ]))
    # select the window
    OClinks_c[min(window.data$order): max(window.data$order),]$window <- i
    
    # stop if there are no more dates to be annotated to a window
    if(length(which(OClinks_c$window == 0)) == 0){break}
    i = max(window.data$order) # start point of next window
  } else {
    next
  }
}

# there are 10 windows of 3 weeks
table(OClinks_c$window)

# number of unique days
days = length(unique(OClinks_c$time))
OClinks_c$day = 0
i = 1; a = 1
for(j in 2:nrow(OClinks_c)){
  day1 = round(OClinks_c$time[i]) 
  day2 = round(OClinks_c$time[j]) 
  if(day1 == day2){
    window.data = OClinks_c[i:j,]
    window.data = unique(rbind(window.data, OClinks_c[OClinks_c$time == tail(window.data, n = 1)$time, ]))
    OClinks_c[min(window.data$order): max(window.data$order),]$day <- a
    i = max(window.data$order) + 1 # start point of next window
    a = a + 1
    if(a > days){break}
  } else {
    next
  }
}

rm(day1); rm(day2); rm(i); rm(j); rm(a); rm(window.data); rm(window)

## for every person, determine:
# average characters send
# which day this person became active (self loop)
# instability of messages to friends over time (is this a person who gets and loses lots of friends?)

attributes$averageChar = 0
attributes$dayactive = 0
attributes$send.day34th.friends = 0
attributes$send.friends = 0

for(i in 1:dim(W)[1]){
  person = OClinks_c[OClinks_c$id1 == i,]
  attributes[i,]$dayactive = person[person$id2 == i,]$day
  attributes[i,]$averageChar = mean(person$characters, na.rm = T)
  
  person.days = length(unique(person$day))
  counter = rep(0, person.days)
  buffer = c()
  for(j in 1:person.days){
    a = length(buffer)
    day.selection = unique(person[person$day == unique(person$day)[j],]$id2)
    buffer = c(buffer, day.selection)
    day.duplicates = length(buffer) - length(unique(buffer))
    buffer = unique(buffer)
    if(day.duplicates != 0){ # not yet in the selection
      counter[j] = length(buffer) - a # add the number of new persons 
    } else { # already yet in the selection
      counter[j] = length(buffer) - a # add the number of new persons 
    }
  }
  
  profile = as.data.frame(cbind(unique(person$day), counter))
  total.contacts = length(unique(person$id2))
  attributes[i,]$send.friends = total.contacts
  profile$sum = cumsum(profile[,2])
  # when does a person reach the 75th percentile or 3/4rds 
  # of his friends
  attributes[i,]$send.day34th.friends = profile[profile$sum >= floor(total.contacts / 4 * 3), 1][1]
  print(i)
}

# small correction if people only send to themselfves (the self-loop to become active)
# the outdegree = 0 so the send.friends should also be zero
attributes[attributes$outdegree == 0, ]$send.friends = 0
attributes[attributes$outdegree == 0, ]$send.day34th.friends = 0

write.table(attributes, "attributes.csv", sep = ";", quote = F, col.names = T, row.names = F)

# Send friends is equal to outdegree in this context so the variable is removed for analyses

# distribution of days when persons reached 75% of their friends
plot(density(attributes[attributes$send.day34th.friends > 0,]$send.day34th.friends),
     main = "Density of days when 75% of total friends had been contacted")

nrow(attributes[attributes$send.day34th.friends > 0,])
plot(attributes[attributes$send.day34th.friends > 0,]$send.day34th.friends,
     1:1350,
     main = "Days when 75% of total friends had been contacted",
     xlab = "time in days",
     ylab = "number of persons")

plot(attributes[attributes$send.day34th.friends > 0,]$send.day34th.friends,
     1:1350,
     main = "Days when 75% of total friends had been contacted",
     type = "l")

# amound and % of persons who had contacted 75% of their friends in the first 50 days
test = attributes[attributes$send.day34th.friends > 0,]
nrow(test[test$send.day34th.friends <= 50,])
nrow(test[test$send.day34th.friends <= 50,]) / nrow(test)

# amound and % of persons who had contacted 75% of their friends in the first 100 days
nrow(test[test$send.day34th.friends <= 100,])
nrow(test[test$send.day34th.friends <= 100,]) / nrow(test)

# amound and % of persons who had contacted 75% of their friends in the first 150 days
nrow(test[test$send.day34th.friends <= 150,])
nrow(test[test$send.day34th.friends <= 150,]) / nrow(test)

# day active
plot(attributes$dayactive, 1:nrow(attributes), xlab = "time in days", ylab = "number of persons")
nrow(attributes[attributes$dayactive <= 50,])
nrow(attributes[attributes$dayactive <= 50,]) / nrow(attributes)

nrow(attributes[attributes$dayactive <= 100,])
nrow(attributes[attributes$dayactive <= 100,]) / nrow(attributes)

nrow(attributes[attributes$dayactive <= 150,])
nrow(attributes[attributes$dayactive <= 150,]) / nrow(attributes)


