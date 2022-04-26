################################################################################
# Dino Dittrich
# 24/04/2017
#
# Program to sample from a posterior based on a normal prior for rho (rho~N(m,std^2)), 
# and a flat prior for beta, and ln(sigma2) for the NAM
# Input in main function (MCMC.N):
#   - # of iterations: R
#   - Vector of observations: y
#   - Design matrix: X
#   - Weight matrix: W
#   - Optional: starting values (as vector for the 3 parameters to be estimated: rho, beta, and sigma2)
# Output:
#   - list containing posterior draws for rho, beta, and sigma2
#   - to be called via "rho.draws","beta.draws",and "sigma2.draws" 
################################################################################
library(mvtnorm)
library(tmvtnorm)
################################################################################
 
# Help functions
# 2 functions for ML computation and MH-step in the MCMC sampler
loglik=function(rho,EW,g,yMy,yMWy,yWMWy){
 (-1)*(sum(log(1-rho*EW))-(g/2)*log(yMy-2*rho*yMWy+rho**2*yWMWy))
}
loglik.gr=function(rho,EW,g,yMy,yMWy,yWMWy){
 (-1)*((-1)*sum(EW/(1-rho*EW))+g*(yMWy-rho*yWMWy)/(yMy-2*rho*yMWy+rho**2*yWMWy))
}
mh.up.rhobeta1.N=function(mu,Sigma,lb,ub,y,Wy,EW,lndet,s2,Xbetatilde,g,SS,cval,Ay,m,std2) {
 pval=c(rtmvnorm(1,mean=mu,sigma=Sigma,lower=c(lb,-Inf),upper=c(ub,Inf)))
 Ay_pval=y-pval[1]*Wy
 lndet_pval=sum(log(1-pval[1]*EW))
 if(log(runif(1))<
  lndet_pval-lndet
  -(1/(2*s2))*(sum(Ay_pval**2)-2*pval[2]*sum(Ay_pval)-2*sum(Ay_pval*Xbetatilde)+g*pval[2]**2+2*pval[2]*sum(Xbetatilde)+sum(Xbetatilde**2)-SS)
  -(1/(2*std2))*(((pval[1]-m)**2)-((cval[1]-m)**2))
  +dmvnorm(cval,mean=mu,sigma=Sigma,log=T)-dmvnorm(pval,mean=mu,sigma=Sigma,log=T)
 ){val=pval;Ay=Ay_pval;lndet=lndet_pval}
 else{val=cval;Ay=Ay;lndet=lndet}
 outlist=list(val,Ay,lndet)
 names(outlist)=c("value","Ay","lndet")
 return(outlist)
}
################################################################################
################################################################################
#
# Main function
#
################################################################################

MCMC.N=function(R=1e4,y,X,W,m,std2,startval){
 g=length(y)
 Id=diag(g)
 ones=rep(1,g)
 Wy=c(W%*%y)
 yWones=sum(Wy)
 yWWy=sum(Wy**2)
 EW=Re(eigen(W)$values)
 tau2=sum(EW**2)
 lb=1/min(EW)
 ub=max(EW)
 k=ncol(X)
 Xtilde=X[,-1]
 XtXi=solve(t(X)%*%X)
 XtXiXt=XtXi%*%t(X)
 M=Id-X%*%XtXiXt
 Sigmai12_help=yWones
 Sigmai21_help=Sigmai12_help
 Sigmai22_help=g
  # Set starting values
  if(missing(startval)) {
  scalepar=.999
  lb_ml=scalepar/min(EW)
  ub_ml=scalepar*max(EW)
  My=M%*%y
  yMy=sum(My**2)
  MWy=M%*%Wy
  yMWy=sum(y*MWy)
  yWMWy=sum(MWy**2)
  rhoS=optim(par=0,fn=loglik,gr=loglik.gr,method="L-BFGS-B",EW=EW,g=g,yMy=yMy,yMWy=yMWy,yWMWy=yWMWy,lower=lb_ml,upper=ub_ml)$par
  Ay=y-rhoS*Wy 
  yAMAy=yMy-2*rhoS*yMWy+rhoS**2*yWMWy
  betaS=c(XtXiXt%*%Ay)
  beta1S=betaS[1]
  betatildeS=betaS[-1]
  s2S=yAMAy/g
  }
  else{
  rhoS=startval[1]
  betaS=startval[2:(1+k)]
  beta1S=betaS[1]
  betatildeS=betaS[-1]
  s2S=startval[2+k]
  Ay=y-rhoS*Wy
  }
 lndet=sum(log(1-rhoS*EW))
 XbetatildeS=c(Xtilde%*%betatildeS)
 SS=sum(Ay**2)-2*beta1S*sum(Ay)-2*sum(Ay*XbetatildeS)+g*beta1S**2+2*beta1S*sum(XbetatildeS)+sum(XbetatildeS**2)
 # Create output vectors and matrices
 rho.draws=NULL
 beta.draws=matrix(NA,nrow=R,ncol=k)
 sigma2.draws=NULL
 rho.draws[1]=rhoS
 beta.draws[1,]=betaS
 sigma2.draws[1]=s2S
  # Run the actual MCMC sampler
  for (r in 2:R) {
   #print(r)
   # Draw (rho,beta1)
   Sigmai11_help=yWWy+s2S*tau2+s2S/std2
   Sigmai_help=rbind(c(Sigmai11_help,Sigmai12_help),c(Sigmai21_help,Sigmai22_help))
   Sigma_help=(1/(Sigmai11_help*Sigmai22_help-Sigmai12_help**2))*rbind(c(Sigmai22_help,-Sigmai12_help),c(-Sigmai21_help,Sigmai11_help))
   Sigma=s2S*Sigma_help
   restilde=c(y-Xtilde%*%betatildeS)
   SStilde=sum(Wy*restilde)
   mu1=sum(restilde*yWones-SStilde-m*s2S/std2)/(yWones**2-g*Sigmai11_help)
    if(mu1>ub|mu1<lb)
    {a=(lb-mu1)/sqrt(Sigma[1,1]);b=(ub-mu1)/sqrt(Sigma[1,1]);z=pnorm(b)-pnorm(a);
    mu1=mu1+(sqrt(Sigma[1,1])/z)*(dnorm(a)-dnorm(b))
    }
   mu2=(SStilde+m*s2S/std2-mu1*Sigmai11_help)/yWones
   mu=c(mu1,mu2)
   Draw=mh.up.rhobeta1.N(mu,Sigma,lb,ub,y,Wy,EW,lndet,s2S,XbetatildeS,g,SS,c(rhoS,beta1S),Ay,m,std2)
   rhoS=Draw$value[1]
   beta1S=Draw$value[2]
   Ay=Draw$Ay
   lndet=Draw$lndet
   # Draw betatilde
   mu_beta=c(XtXiXt%*%Ay)
   Sigma_beta=s2S*XtXi
   Sigma_beta11=c(Sigma_beta[1:1])
   Sigma_beta12=c(Sigma_beta[1,2:k])
   Sigma_beta22=Sigma_beta[2:k,2:k]
   mu_betatilde=mu_beta[-1]+Sigma_beta12*(beta1S-mu_beta[1])/Sigma_beta11
   Sigma_betatilde=Sigma_beta22-outer(Sigma_beta12,Sigma_beta12)/Sigma_beta11
   betatildeS=c(rmvnorm(1,mean=mu_betatilde,sigma=Sigma_betatilde))
   XbetatildeS=c(Xtilde%*%betatildeS)
   SS=sum(Ay**2)-2*beta1S*sum(Ay)-2*sum(Ay*XbetatildeS)+g*beta1S**2+2*beta1S*sum(XbetatildeS)+sum(XbetatildeS**2)
   # Draw sigma2
   s2S=1/rgamma(1,shape=g/2,rate=SS/2)
   rho.draws[r]=rhoS
   beta.draws[r,]=c(beta1S,betatildeS)
   sigma2.draws[r]=s2S
  }
 # Return the posterior draws
 outlist=list(rho.draws,beta.draws,sigma2.draws)
 names(outlist)=c("rho.draws","beta.draws","sigma2.draws")
 return(outlist)
}
################################################################################
################################################################################
#
# Example for n=50
#
################################################################################
library(sna)
library(numDeriv)
################################################################################

# set network size (n), network density(d), and generate random W
set.seed(2)
n=50
d=.05
Wadj=rgraph(n,tprob=d,mode="graph")
W=make.stochastic(Wadj,mode="row")
# set rho, beta, sigma2, and generate y
rho=.05
beta=rnorm(4)
sigma2=1
X=matrix(c(rep(1,n),rnorm(n*3)),n,4)
y=c((solve(diag(n)-rho*W))%*%(X%*%beta+rnorm(n)))
# set prior mean and standard deviation (for rho)
m=.36; std2=.19**2
# run MCMC sampler
out=MCMC.N(R=1e4,y,X,W,m,std2)
# run MLE
out.mle=lnam(y,X,W)
# posterior median estimates
c(median(out$rho.draws),apply(out$beta.draws,2,median),median(out$sigma2.draws)) 
# ML estimates
c(out.mle$rho1,out.mle$beta,out.mle$sigma**2)

#####
# trace plot of posterior draws for rho with data-generating value (red) and prior mean (green) added
plot(out$rho.draws,type="l",ylab=expression(rho))
abline(h=rho,col=2)
abline(h=m,col=3)
# empirical autocorrelation function for posterior draws for rho
acf(out$rho.draws)
# marginal posterior density of rho with data-generating value (red) and prior mean (green) added
plot(density(out$rho.draws),xlab=expression(rho),main="")
abline(v=rho,col=2)
abline(v=m,col=3)
# marginal posterior density of first regression coefficient (not the intercept) with data-generating value added
plot(density(out$beta.draws[,2]),xlab=expression(beta[1]),main="")
abline(v=beta[2],col=2)
# marginal posterior density of sigma2 with data-generating value added
plot(density(out$sigma2.draws),xlab=expression(sigma^2),main="")
abline(v=sigma2,col=2)



