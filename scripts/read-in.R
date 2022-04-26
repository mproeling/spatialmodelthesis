rm(list = ls(all = TRUE))
set.seed(2017)
library(spatialprobit)

data("Katrina")
attach(Katrina) 
table(y1) # 300 of the 673 firms reopened during 0-3 months horizon, p.1016 
table(y2) # 425 of the 673 firms reopened during 0-6 months horizon, p.1016 
table(y3) # 478 of the 673 firms reopened during 0-12 months horizon, p.1016 detach(Katrina)
detach(Katrina)

if (require(ggmap)) { qmplot(long, lat, data = Katrina, maptype="roadmap", source="google") }
require(spdep)

# (a) 0-3 months time horizon # LeSage et al. (2011) use k=11 nearest neighbors in this case 
nb <- knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k=11)) 
listw <- nb2listw(nb, style="W")
W1 <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
detach("package:spatialprobit", unload = TRUE)

#DataSet = as.data.frame(cbind(1,Katrina[,c(4,5,6,7,8,12)]))
#colnames(DataSet)[1] = c("Intercept")
X = as.matrix(cbind(1, Katrina[,c(4,5,6,7,8,9,10,11)]))
colnames(X)[1] = "Intercept"
#y = rnorm(nrow(Katrina))
y = Katrina[,12]
#y = cbind(Katrina[,12], Katrina[,13]) 
#y = cbind(y, rnorm(nrow(y), 0, 1))

model1 = sarprobit(y1 ~ flood_depth + small_size + large_size, W = W1, ndraw=1000, burn.in = 100, showProgress = TRUE, data = Katrina)

ndraw=10000
burn.in=1000
thinning=1
W = W1
m = 10
rho = 0.39
rhosd = 0.16
x=X

beta = list(); beta[[1]] = rep(0, ncol(X)) # for 1 y variable
#beta = list(); for(i in 1:ncol(y)){beta[[i]] = rep(0, (ncol(X) + (i-1)))}

computeMarginalEffects=TRUE
showProgress=TRUE

source('/Users/mark/Documents/networks/scripts/gibbs_zsamplingfunctions.R')
source('/Users/mark/Documents/networks/scripts/sar_base.r')
source('/Users/mark/Documents/networks/scripts/sarprobit_sarcontinuous_functions.R')

prior2 = c()
prior2$novi = 1

# create missing data block
# for 1 Y covariate
y[23]<-NA # create missingness where Y = 0
#y[381]<-NA # create missingness where Y = 1
#X[c(23, 381),] 
method = c("continuous")

# for 2 Y covariates
#y[23,] <- NA
#y[381,] <- NA
#y[45,2] <- NA
#y[250,2] <- NA
#method = c("probit", "continuous")

# for 3 Y covariates
y[23,] <- NA
y[381,] <- NA
y[45,c(2,3)] <- NA
y[250,c(2,3)] <- NA
y[115,c(3)] <- NA
y[120,c(3)] <- NA
method = c("probit", "probit", "continuous")

prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1)

#model = sar_combined_mcmc_imp_cut(y, x, W, ndraw=10000, burn.in=1000, thinning=1, 
#                      prior=prior, 
#                      rho = 0.39, rhosd = 0.16, 
#                      m=10, computeMarginalEffects=TRUE, showProgress=TRUE,
#                      method = method)



#######################################################
#######################################################
#######################################################

rm(list = ls(all = TRUE))
set.seed(2017)
library(spatialprobit)

load("C:/Users/mark/Downloads/Glasgow_data/Glasgow-friendship.RData")
load("C:/Users/mark/Downloads/Glasgow_data/Glasgow-demographic.RData")
load("C:/Users/mark/Downloads/Glasgow_data/Glasgow-geographic.RData")
load("C:/Users/mark/Downloads/Glasgow_data/Glasgow-lifestyle.RData")
load("C:/Users/mark/Downloads/Glasgow_data/Glasgow-selections.RData")
load("C:/Users/mark/Downloads/Glasgow_data/Glasgow-substances.RData")
load("C:/Users/mark/Downloads/Glasgow_data/Glasgow-various.RData")

# remove the persons with value 10 in adjacency matrix
include = which(friendship.1[,1] != 10)
friendship.1.clean = friendship.1[include,include]
table(friendship.1.clean==10)
table(friendship.1.clean)
# recode friendschip weights:
# 0 = no friends -> 0
# 1 = best friends -> 1 
# 2 = just friends -> 0.5

for(i in 1:nrow(friendship.1.clean)){
  for(j in 1:ncol(friendship.1.clean)){
   if(friendship.1.clean[i,j] == 2){friendship.1.clean[i,j] = 0.5}
  }
}
table(friendship.1.clean)

alcohol.1 = alcohol[include,]
cannabis.1 = cannabis[include,]
tobacco.1 = tobacco[include,]
familysmoking.1 = familysmoking[include,]
money.1 = money[include,]
leisure.1 = leisure1[include,]
music.1 = music1[include,]
romantic.1 = romantic[include,]
age.1 = age[include]
sex.1 = sex.F[include]









