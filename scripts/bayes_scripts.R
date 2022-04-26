rm(list = ls(all = TRUE))
library(MASS)

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lmfit <- lm(weight ~ group)

## function to compute the bayesian analog of the lmfit
## using non-informative priors and Monte Carlo scheme
## based on N samples

N = 1e4

QR <- lmfit$qr
df.residual <- lmfit$df.residual
R <- qr.R(QR) ## R component
coef <- lmfit$coef
Vb <- chol2inv(R) ## variance(unscaled)
s2 <- (t(lmfit$residuals)%*%lmfit$residuals)
s2 <- s2[1,1]/df.residual

## now to sample residual variance
sigma <- df.residual*s2/rchisq(N,df.residual)
coef.sim <- sapply(sigma,function(x) mvrnorm(1,coef,Vb*x))
ret <- data.frame(t(coef.sim))
names(ret) <- names(lmfit$coef)
ret$sigma <- sqrt(sigma)
ret

hist(coef.sim[1,])
hist(coef.sim[2,])
hist(ret[,3])

Bayes.sum<-function(x)
{
  c("mean"=mean(x),
    "se"=sd(x),
    "t"=mean(x)/sd(x),
    "median"=median(x),
    "CrI"=quantile(x,prob=0.025),
    "CrI"=quantile(x,prob=0.975)
  )
}

set.seed(1234)  ## reproducible sim
lmfit <- lm(weight ~ group)
bf<-bayesfit(lmfit,10000)
t(apply(bf,2,Bayes.sum))
summary(lmfit)

######
a= cbind(c(90,90,60,30,30),
        c(80,60,50,40,20),
        c(40,80,70,70,90))
b = matrix(rep(1,25), 5, 5)

A = a - b%*%a*1/5
ATA = t(A)%*%A




A= cbind(c(4,12,-16),
         c(12,37,-43),
         c(-16,-43,98))

ca = chol(A)
chol2inv(ca)

#### optim example
dat=data.frame(x=c(1,2,3,4,5,6), 
               y=c(1,3,5,6,8,12))

min.RSS <- function(data, par) {
  with(data, sum((par[1] + par[2] * x - y)^2))
}

model = optim(par = c(0,1), min.RSS, data = dat)

plot(y ~ x, data = dat)
abline(a = model$par[1], b = model$par[2], col = "red")

###### logistic bayes
"There are several packages in R that include MCMC approaches. Here we use the MCMCpack package, which include the MCMClogit() function.
It appears not to accept the weights option mentioned previously, so we generate data at the observation level to begin. Then we run the MCMC."
events.0=0   # for X = 0
events.1=5   # for X = 1
x = c(rep(0,100), rep(1,100))
y = c(rep(0,100-events.0), rep(1,events.0),
      rep(0, 100-events.1), rep(1, events.1))

library(MCMCpack)
logmcmc = MCMClogit(y~as.factor(x), burnin=1000, mcmc=21000, b0=0, B0=.04)

"The MCMClogit() accepts a formula object and allows the burn-in and number of Monte Carlo iterations to be specified. 
The prior mean b0 can be specified as a vector if different for each parameter, or as a scalar, as show. Similarly,
the prior precision B0 can be a matrix or a scalar; if scalar, the parameters are uncorrelated in the prior."

summary(logmcmc)
plot(logmcmc)

##################################################################################
#   from: 
#   https://github.com/stablemarkets/BayesianTutorials/tree/master/SimpleLinearReg
##################################################################################

rm(list = ls(all = TRUE))

library(invgamma)
library(MASS)
library(xtable)
set.seed(1)

#############################################################################
###########   0 - Simulate Data 
#############################################################################
# beta_0 hyper parameters (known)
m0<-5
t0<-1

m1<-5
t1<-1

a<-.5 # shape
g<-.7 # scale

n<-100
x<-rnorm(n, 0, 1)

tphi<-rinvgamma(1, shape=a, rate=g)
tb0<-rnorm(1, m0, sqrt(t0) ) # intercept
tb1<-rnorm(1, m1, sqrt(t1) ) # betas
tphi; tb0; tb1;

y<-rnorm(n, tb0 + tb1*x, sqrt(tphi))

#############################################################################
###########   1 - Functions for grid evaluation of posterior densities
#############################################################################
# b0 conditional
rb0cond<-function(y, x, b1, phi, t0, m0){
  grid<-seq(-10,10,.001)
  
  p<-numeric(length=length(grid))
  for(i in 1:length(p) ){
    p[i]<- (-(1/(2*phi))*sum( (y - (grid[i]+b1*x))^2 ))  + ( -(1/(2*t0))*(grid[i] - m0)^2)
  }
  
  draw<-sample(grid, size = 1, prob = exp(1-p/max(p)))
  return(draw)
}

# b1 conditional
rb1cond<-function(y, x, phi, t1, m1, b0){
  grid<-seq(-10,10,.001)
  
  p<-numeric(length=length(grid))
  for(i in 1:length(p)){
    p[i]<- (-(1/(2*phi) )*sum( (y - (b0+grid[i]*x))^2 )) + ( -(1/(2*t1) )*(grid[i] - m1)^2)
  }
  
  draw<-sample(grid, size = 1, prob = exp(1-p/max(p)))
  return(draw)
}

#############################################################################
###########   2 -Implement Gibbs Sampling
#############################################################################

iter<-1000
burnin<-101
phi<-b0<-b1<-numeric(iter)
phi[1]<-b0[1]<-b1[1]<-6

for(i in 2:iter ){
  # taking as pho the previous run of b0 and b1
  phi[i]<-rinvgamma(1, shape = (n/2 + a), rate = .5*sum( (y - (b0[i-1]+b1[i-1]*x))^2 ) + g  )
  b0[i]<-rb0cond(y, x, b1[i-1], phi[i], t0, m0) # taking the previous run of b1
  b1[i]<-rb1cond(y, x, phi[i], t1, m1, b0[i] )
}

#############################################################################
###########   3 - Visualize Results
#############################################################################

par(mfrow=c(2,2))
plot(phi[burnin:iter],type='l'); abline(h=tphi, col='red')
plot(b0[burnin:iter],type='l'); abline(h=tb0, col='red')
plot(b1[burnin:iter],type='l'); abline(h=tb1, col='red')

z <- kde2d(b0, b1, n=50)
plot(b0,b1, pch=19, cex=.4)
contour(z, drawlabels=FALSE, nlevels=10, col='red', add=TRUE)
