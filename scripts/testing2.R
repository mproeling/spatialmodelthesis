## Not run:
rm(list = ls(all = TRUE))

library(spatialprobit)
set.seed(1)
################################################################################
#
# Example with J = 4 alternatives
#
################################################################################
# set up a model like in SAR probit
J <- 4
# ordered alternatives j=1, 2, 3, 4
# --> (J-2)=2 cutoff-points to be estimated phi_2, phi_3
phi <- c(-Inf, 0, +1, +2, Inf) # phi_0,...,phi_j, vector of length (J+1)
# phi_1 = 0 is a identification restriction

# generate random samples from true model
n <- 400 # number of items
k <- 3 # 3 beta parameters
beta <- c(0, -1, 1) # true model parameters k=3 beta=(beta1,beta2,beta3)
rho <- 0.75
# design matrix with two standard normal variates as "coordinates"
X <- cbind(intercept=1, x=rnorm(n), y=rnorm(n))
# identity matrix I_n
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
# build spatial weight matrix W from coordinates in X
W <- kNearestNeighbors(x=rnorm(n), y=rnorm(n), k=6)
# create samples from epsilon using independence of distributions (rnorm())
# to avoid dense matrix I_n
eps <- rnorm(n=n, mean=0, sd=1)
z <- solve(qr(I_n - rho * W), X)
# ordered variable y:
# y_i = 1 for phi_0 < z <= phi_1; -Inf < z <= 0
# y_i = 2 for phi_1 < z <= phi_2
# y_i = 3 for phi_2 < z <= phi_3
# y_i = 4 for phi_3 < z <= phi_4# y in {1, 2, 3}
y <- cut(as.double(z), breaks=phi, labels=FALSE, ordered_result = TRUE)
table(y)
#y
# 1 2 3 4
#246 55 44 55
# estimate SAR Ordered Probit
res <- sar_ordered_probit_mcmc(y=y, X=X, W=W, showProgress=TRUE)
summary(res)
           