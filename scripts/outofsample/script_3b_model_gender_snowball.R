rm(list = ls(all = TRUE))

#missing_person_id = 1 
#args=(commandArgs(TRUE))
#missing_person_id = eval(parse(text=args[[1]]))
#missing_person_id = 1 

# this script is updated by taking Year of Study as a set of dummy variables
library(igraph)
library(spatialprobit)
library(corrplot)
library(truncnorm)
library(snowboot)

source('D:/Documents/networks/scripts/outofsample/gibbs_zsamplingfunctions.R')
source('D:/Documents/networks/scripts/outofsample/sar_base.R')
source('D:/Documents/networks/scripts/outofsample/sar_continuous.R')

function.list = ls(all = TRUE)

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

beta = list(); beta[[1]] = rep(0, ncol(X)) # for 1 y variable
#beta = rep(0, ncol(X))
computeMarginalEffects=TRUE
showProgress=TRUE

prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1)
start=list(rho=0.39, rhosd=0.16, beta=rep(0, ncol(X)), phi=c(-Inf, 0:(max(y)-1), Inf))

prior2 = c()
prior2$novi = 1

##################################################
#################################################
library(parallel)
library(foreach)
library(doParallel)

source("D:/Documents/networks/scripts/outofsample/sar_imp_cut_oneout.R")
missing.functie <- function(x,y){
  yf1 = y
  yf1[unique(x)] <- NA
  sar_combined_mcmc_imp_cut(y = yf1, x=X, W, ndraw=1000, burn.in=50, thinning=1, start = start,
                            prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                            m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                            model = x)
}

setwd("D:\\Documents\\downloads\\ToreOpsahl\\output")

cores = detectCores()
cores = 7 
clust <- makeCluster(cores)
doParallel::registerDoParallel(clust)

clusterExport(clust, list("missing.functie", "X", "y", "W", "start", "prior", "prior2", "beta", "sar_combined_mcmc_imp_cut"))
clusterExport(clust, function.list)
clusterEvalQ(clust, library(igraph))
clusterEvalQ(clust, library(spatialprobit))
clusterEvalQ(clust, library(corrplot))
clusterEvalQ(clust, library(truncnorm))
clusterEvalQ(clust, library(snowboot))

#output = parLapply(clust, 1:5, function(x) missing.functie(x, y))

for(i in 1:5){
  a = ((i-1) * 50) + 1
  b = i * 50
  output = parLapply(clust, a:b, function(x) missing.functie(x, y))
  print(i)
}

stopCluster(clust)
          






          
#### cut model
source("D:/Documents/networks/scripts/outofsample/sar_imp_cut_oneout.R")
Sys.time()
model = sar_combined_mcmc_imp_cut(y = yf1, x=X, W, ndraw=1000, burn.in=100, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = x)
Sys.time()
#write.table(impmodel.cluster.normal.10$y_imp[[1]], "y_imp_normalcut.csv", sep = ";", quote = F, col.names = F, row.names = F)

source("D:/Documents/networks/scripts/outofsample/sar_imp_fbs_oneout.R")
model = sar_combined_mcmc_imp_fbs(y = yf1, x=X, W, ndraw=1000, burn.in=100, thinning=1, start = start,
                                        prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                        m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                        model = missing_person_id)
#write.table(impmodel.cluster.normal.fbs.10$y_imp[[1]], "y_imp_normalfbs.csv", sep = ";", quote = F, col.names = F, row.names = F)



