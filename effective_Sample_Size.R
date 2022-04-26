library(coda)
library(mcmcse)

setwd("D:/Documents/downloads/ToreOpsahl/bayes_oxford_snowball_gender_18052018/rhoestimates/")
setwd("D:/Documents/downloads/ToreOpsahl/bayes_snowball_gender_28052018")
# read in the rho estimates for the full data data
effective0 = read.table("rhoestimatesnowball/rhoestimate_nomissing_cut.txt")

# read in the rho estimates with missing data
effective10cut = read.table("rhoestimatesnowball/rhoestimate10cut.txt")
effective10fbs = read.table("rhoestimatesnowball/rhoestimate10fbs.txt")

effective25cut = read.table("rhoestimatesnowball/rhoestimate25cut.txt")
effective25fbs = read.table("rhoestimatesnowball/rhoestimate25fbs.txt")

effective50cut = read.table("rhoestimatesnowball/rhoestimate50cut.txt")
effective50fbs = read.table("rhoestimatesnowball/rhoestimate50fbs.txt")

effective75cut = read.table("rhoestimatesnowball/rhoestimate75cut.txt")
effective75fbs = read.table("rhoestimatesnowball/rhoestimate75fbs.txt")

test = as.mcmc(effective0$V1)
test1 = as.mcmc(effective10cut$V1)

lapply(test, effectiveSize)



function (x) 
{
  if (is.mcmc.list(x)) {
    ess <- do.call("rbind", lapply(x, effectiveSize))
    ans <- apply(ess, 2, sum)
  }
  else {
    x <- as.mcmc(x)
    x <- as.matrix(x)
    spec <- spectrum0.ar(x)$spec
    ans <- ifelse(spec == 0, 0, nrow(x) * apply(x, 2, var)/spec)
  }
  return(ans)
}



function (x) 
{
  x <- as.matrix(x)
  v0 <- order <- numeric(ncol(x))
  names(v0) <- names(order) <- colnames(x)
  z <- 1:nrow(x)
  for (i in 1:ncol(x)) {
    lm.out <- lm(x[, i] ~ z)
    if (identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
      v0[i] <- 0
      order[i] <- 0
    }
    else {
      ar.out <- ar(x[, i], aic = TRUE)
      v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
      order[i] <- ar.out$order
    }
  }
  return(list(spec = v0, order = order))
}



#x <- c(rnorm(1e3), rnorm(1e3, 10))
#ess(x1)
#x1 <- x[1:1000]
#x2 <- x[1001:2000]
#multiESS(cbind(x1,x2))

x0 <- effective0[c(1001:26000),1]
x1 <- effective10cut[c(1001:26000),1]
x2 <- effective25cut[c(1001:26000),1]
x3 <- effective50cut[c(1001:26000),1]
x4 <- effective75cut[c(1001:26000),1]

x5 <- effective10fbs[c(1001:26000),1]
x6 <- effective25fbs[c(1001:26000),1]
x7 <- effective50fbs[c(1001:26000),1]
x8 <- effective75fbs[c(1001:26000),1]


multiESS(cbind(x0,x2))
multiESS(cbind(x0,x3))
multiESS(cbind(x0,x4))
multiESS(cbind(x0,x5))
multiESS(cbind(x0,x6))
multiESS(cbind(x0,x7))
multiESS(cbind(x0,x8))

# voor 75%
x0 <- effective0[c(1001:26000),1]
x1 <- effective75cut[c(1001:26000),1]
x2 <- effective75fbs[c(1001:26000),1]

multiESS(cbind(x0,x1))
multiESS(cbind(x0,x2))

ess(x1)
ess(x2)
ess(x3)
ess(x4)
ess(x5)
ess(x6)
ess(x7)
ess(x8)

multiESS(cbind(x1,x5))
multiESS(cbind(x2,x6))
multiESS(cbind(x3,x7))
multiESS(cbind(x4,x8))


