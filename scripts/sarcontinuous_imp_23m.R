###########################################################
##### Bayesian Spatial regression model 
##### Continuous imputation
###########################################################
rm(list = ls(all = TRUE))

source('/Users/mark/Documents/networks/scripts/functions/gibbs_zsamplingfunctions.R')
source('/Users/mark/Documents/networks/scripts/functions/sar_base.r')
library(igraph)
set.seed(2017)

anselin = read.table("/Users/mark/Downloads/jplv7/data/anselin.dat")
colnames(anselin) = c("crime", "income", "values", "latt", "long")
n = nrow(anselin)
y = anselin[,1] # select crime

X = cbind(rep(1, n), anselin[,2:3]) 
latt = anselin[,4]
long = anselin[,5]

# 1st order contiguity matrix for
# Anselin's Columbus crime dataset
# stored in sparse matrix format [i, j, s] = find(W);
# so that W = sparse(i,j,s); reconstructs the 49x49 matrix
# NOTE: already row-standardized

W = read.table("/Users/mark/Downloads/jplv7/data/wmat.dat")
W[,1]=as.character(W[,1]) #Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
W[,2]=as.character(W[,2])
W=as.matrix(W) #igraph needs the edgelist to be in matrix format
g=graph.edgelist(W[,1:2]) #We first greate a network from the first two columns, which has the list of vertices
E(g)$weight=as.numeric(W[,3]) #We then add the edge weights to this network by assigning an edge attribute called 'weight'.
W = as_adjacency_matrix(g, attr="weight", sparse = TRUE)
W = W[order(as.numeric(colnames(W))),order(as.numeric(colnames(W)))]

info.lflag = 0

prior2 = c()
prior2$novi = 1 # homoscedastic model
prior2$lflag = 0

sar_combined_mcmc_imp_fbs = function(y, x, W, ndraw=1000, burn.in=100, thinning=1, 
                               prior=list(a1=1, a2=1, beta = as.list(beta), T=diag(ncol(X))*1e12, lflag = 0), 
                               rho = 0.39, rhosd = 0.16, 
                               m=10, computeMarginalEffects=TRUE, showProgress=FALSE,
                               method = vector("character", meth.length = ncol(x))){}

# definitie van nieuwe imputatie functie met echte data

  
  
  



