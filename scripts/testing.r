rm(list = ls(all = TRUE))
source('/Users/mark/Documents/networks/scripts/gibbs_zsamplingfunctions.R')
source('/Users/mark/Documents/networks/scripts/sar_base.r')
source('/Users/mark/Documents/networks/scripts/sarprobit_sarcontinuous_functions.R')
library(igraph)
library(spatialprobit)

data("Katrina")
attach(Katrina) 
table(y1) # 300 of the 673 firms reopened during 0-3 months horizon, p.1016 
table(y2) # 425 of the 673 firms reopened during 0-6 months horizon, p.1016 
table(y3) # 478 of the 673 firms reopened during 0-12 months horizon, p.1016 detach(Katrina)
detach(Katrina)
# detach("package:spatialprobit", unload=TRUE)

if (require(ggmap)) { qmplot(long, lat, data = Katrina, maptype="roadmap", source="google") }
require(spdep)
# (a) 0-3 months time horizon # LeSage et al. (2011) use k=11 nearest neighbors in this case 
nb <- knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k=11)) 
listw <- nb2listw(nb, style="W")
W1 <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")

fit1 <- sarprobit(y1 ~ flood_depth + log_medinc + small_size + large_size + low_status_customers + 
                    high_status_customers + owntype_sole_proprietor + owntype_national_chain, 
                  W=W1, data=Katrina, ndraw=600, burn.in = 100, showProgress=TRUE)
summary(fit1)

###################################### 
##### TEST CONTINOUS FUNCTION  
###################################### 

rm(list = ls(all = TRUE))
source('/Users/mark/Documents/networks/scripts/functions/gibbs_zsamplingfunctions.R')
source('/Users/mark/Documents/networks/scripts/functions/sar_base.r')
library(igraph)

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


model1 = sar_continuous_mcmc(y, X, W, ndraw=2500, burn.in=500, thinning=1, 
                             prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
                             m=10, computeMarginalEffects=TRUE, showProgress=TRUE)

# compare my function against the script from Dino Dittrich (adapted to accomodate the sparse W matrix)
source('/Users/mark/Documents/networks/scripts/NAM_MCMC_continuous_mp.R')
model1.Dino = out=MCMC.N(R=1e4,y,X,W,m,std2)


###############################################################
# continue script by fitting model with heteroscedastic prior
###############################################################

model1$tflag = 'tstat'

# Gibbs sampling function heteroscedastic prior
# to maximum likelihood estimates
prior = c()
prior$rval = 4
prior2$novi = 0 # heteroscedastic model
prior$lflag = 0

model2 = sar_continuous_mcmc(y, X, W, ndraw=2500, burn.in=500, thinning=1, 
                             prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
                             m=10, computeMarginalEffects=TRUE, showProgress=FALSE)

###############################################################
# run imputation model
###############################################################
info.lflag = 0

prior2 = c()
prior2$novi = 1 # homoscedastic model
prior2$lflag = 0

y[23]<-NA # create missingness

# model with beta priors defined
model3a = sar_continuous_mcmc_imp(y, X, W, ndraw=2500, burn.in=500, thinning=1, 
                                 prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
                                 m=10, computeMarginalEffects=TRUE, showProgress=TRUE)

# model with beta prior not defined
model3b = sar_continuous_mcmc_imp(y, X, W, ndraw=2500, burn.in=500, thinning=1, 
                                 prior=list(a1=1, a2=1, T=diag(ncol(X))*1e12, lflag = 0), 
                                 m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("pussy","pussy"))


# multiple y with missingness
y2 = y 
y2 = rnorm(n, y2, sd(y2))
y = cbind(y,y2)
y[23,]<-NA # create missingness
y[c(10,11,12),2] <- NA # create missingness

model4 = sar_continuous_mcmc_imp(y, X, W, ndraw=2500, burn.in=500, thinning=1, 
                                 prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
                                 m=10, computeMarginalEffects=TRUE, showProgress=FALSE)

# comparison with normal linear model
lm(y2 ~ y[,1] + X[,2] + X[,3])


# second test
y[,2] = y[,1]
y[10,2] = NA
model5 = sar_continuous_mcmc_imp(y, X, W, ndraw=2500, burn.in=500, thinning=1, 
                                 prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0), 
                                 m=10, computeMarginalEffects=TRUE, showProgress=FALSE)
colMeans(model5$y_imp[[1]])
colMeans(model5$y_imp[[2]])

