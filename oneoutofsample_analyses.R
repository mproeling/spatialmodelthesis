rm(list = ls(all = TRUE))

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

# load imputed data for every observation
folder1 = "D:/Documents/downloads/ToreOpsahl/outofsample_cut/outofsample/"
folder2 = "D:/Documents/downloads/ToreOpsahl/outofsample_fbs/outofsample/"
cutpath = folder1
fbspath = folder2
output = matrix(0, 1899, 2)
file.names.cut <- dir(cutpath, pattern ="imp_")
file.names.fbs <- dir(fbspath, pattern ="imp_")

for(i in 1:length(y)){

  ## GET DATA
  #get y value
  y.value = y[i]
  
  #cutmodel
  cutfile <- read.table(paste0(folder1, file.names.cut[i]))
  cutfile = as.data.frame(cutfile[-c(1:50),]) # drop burn in period
  colnames(cutfile) = "V1"
  
  #fbsmodel
  fbsfile <- read.table(paste0(folder2, file.names.fbs[i]))
  fbsfile = as.data.frame(fbsfile[-c(1:50),]) # drop burn in period
  colnames(fbsfile) = "V1"
  
  ## ANALYSE
  cutfile$sign = ifelse(cutfile$V1 <= 0, 0, 1)
  fbsfile$sign = ifelse(fbsfile$V1 <= 0, 0, 1)
  
  cutmodelwins = NA
  if(y.value == 1){
    cutmodelwins = ifelse(sum(cutfile$sign) > sum(fbsfile$sign), 1, 0)
  } else {
    cutmodelwins = ifelse(sum(cutfile$sign) < sum(fbsfile$sign), 1, 0)
  }
  
  output[i,1] = i
  output[i,2] = cutmodelwins
  print(i)
}


