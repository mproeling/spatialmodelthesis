# donor selection function

donor_selection = function(W, y, y_pred, imp.obs, Nd = 5, Nmiss, missing){

  donor.output = matrix(0, Nd, 1) 
  for(MIS in 1:Nmiss[[i]]){
    person_i = missing[[i]][MIS] # id of person i
    W_i = cbind(seq(1, nrow(W)), W[,person_i]) # get W of person i
    n1i = sort(W_i[W_i[,2] != 0, 1]) # # get W of person i
    # select \hat(y) data from first degree neighbours of person i
    x_n1i = y_pred[n1i,]
    
    # calculate difference of \hat(y) and \hat(y)_{obs}
    imp.obs$diff = abs(imp[MIS,]$y_pred - imp.obs$y_pred)
    imp.obs = imp.obs[order(imp.obs$diff),]
    
    output = data.frame()
    # for every person i, select a donor
    for(u_d in 1:Nd){
      # get donor id from sorted imp.obs 
      donor = sort(imp.obs[u_d,]$id)
      # select the neighbour of donors, get W
      Wdonor = cbind(seq(1, nrow(W)), matrix(W[,donor]))
      n1u = sort(unique(Wdonor[Wdonor[,2] != 0, 1]))
      x_n1u = y_pred[n1u,]
      x_n1iu = rbind(x_n1i, x_n1u)
      
      person.data = y_pred[y_pred$id == person_i,]
      donor.data = y_pred[y_pred$id == donor,]
      
      output[u_d,1] = person_i      # id person i
      output[u_d,2] = nrow(x_n1i)   # neighbours person i
      output[u_d,3] = donor         # id donor
      output[u_d,4] = nrow(x_n1u)   # neighbours donor
      
      # 1a: find the shared nodes
      n_occur <- data.frame(table(x_n1iu$id))
      shared = as.data.frame(n_occur[which(n_occur$Freq > 1), 1])
      not.shared = as.data.frame(n_occur[which(n_occur$Freq < 2), 1])
      colnames(not.shared) = "id"
      colnames(shared) = "id"
      if(nrow(shared) > 0){
        output[u_d,5] = length(n_occur[n_occur$Freq > 1,]$Var1) # amount of shared nodes
        
        # 1b: ratio between d(shared nodes, i) / d(shared nodes, u)
        x_n1iu = unique(rbind(x_n1i, x_n1u))
        shared.data = merge(x_n1iu, shared)
        combined.data = cbind(shared.data, person.data, donor.data)
        combined.data$person.diff = abs(combined.data[,2] - combined.data[,4]) # difference in y_pred shared friend(s) - person
        combined.data$donor.diff = abs(combined.data[,2] - combined.data[,6]) # difference in y_pred shared friend(s) - donor
        output[u_d,6] = mean(combined.data$person.diff) / mean(combined.data$donor.diff)
        
      } else {
        output[u_d,5] = 0
        output[u_d,6] = NA
      }
      
      # 2: d(i, nodes in cluster i) / d(u, nodes in cluster u), but not the similar nodes 
      donor.friends =  merge(x_n1u, not.shared)
      donor.friends = donor.friends$y_pred
      if(length(donor.friends) == 0){
        output[u_d,7] = NA
        output[u_d,8] = NA
        next
      } else {
        donor.cluster = as.data.frame(cbind(donor.data[,2], donor.friends))
        donor.cluster$diff = abs(donor.cluster[,1] - donor.cluster[,2])
        person.friends = merge(x_n1i, not.shared)
        person.friends = person.friends$y_pred
      }
      
      if(length(person.friends) == 0){
        output[u_d,7] = NA
        output[u_d,8] = NA
        next 
      } else {
        person.cluster = as.data.frame(cbind(person.data[,2], person.friends))
        person.cluster$diff = abs(person.cluster[,1] - person.cluster[,2])
        output[u_d,7] = mean(person.cluster$diff) / mean(donor.cluster$diff)
      }
      
      # 3: d(shared nodes van i - shared nodes van u)
      distance = matrix(0, length(person.friends) * length(donor.friends), 1)
      count = 1 
      for(ab in 1:length(person.friends)){
        for(ba in 1:length(donor.friends)){
          distance[count,] = abs(person.friends[ab] - donor.friends[ba])
          count = count + 1 
        }
      }
      output[u_d,8] = mean(distance)
    }
    
    # hoe meer shared friends hoe beter
    output_Nneighbours = log(output[,c(4,5)])
    output_Nneighbours[output_Nneighbours == "-Inf"] <- 0
    output_Nneighbours_weights = rowSums(output_Nneighbours)
    output_Nneighbours_weights[output_Nneighbours_weights == 0] <- 1
    
    # 6 = verhouding mean(difference(shared friends - person)) / mean(difference(shared friends - donor))
    # 7 = verhouding mean(difference(non shared friends - person)) / mean(difference(non shared friends - donor))
    output_difference_ratio_1 = rank(1 - output[,6])
    output_difference_ratio_2 = rank(1 - output[,7])
    
    # 8 is the mean of the difference between all persons of the 2 two clusters 
    output_difference = rank(output[,8])
    
    donor_weights = (output_difference_ratio_1 + output_difference_ratio_2 + output_difference) / output_Nneighbours_weights
    # select donor 
    donor_id = output[which(donor_weights == min(donor_weights)), 3]
    donor_values = imp.obs[imp.obs$id == donor_id,]$y_pred

    # output = cbind(output, donor_id)
    # we now finished donor - friends - person distance calculation and can select a donor 
    # first write the output
    if(iter == 1 - burn.in){
      write.table(output, paste0("clusteroutput", outputname, ".txt"), col.names = F, row.names = F, quote = F, sep = "\t")
    } else {
      write.table(output, paste0("clusteroutput", outputname, ".txt"), col.names = F, row.names = F, quote = F, sep = "\t", append = T)
    }
    
    donor.output[Nd,] = mean(donor_values)
  }
  return(donor.output)
}

  