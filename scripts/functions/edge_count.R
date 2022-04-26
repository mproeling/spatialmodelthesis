
sar_combined_mcmc_imp_cut = function(y, x, W, ndraw=1000, burn.in=100, thinning=1, 
                                     prior=list(a1=1, a2=1, beta = as.list(beta), T=diag(ncol(X))*1e12, lflag=0), 
                                     start=list(rho=0.39, rhosd=0.16, beta=rep(0, ncol(X)), phi=c(-Inf, 0:(max(y)-1), Inf)),
                                     m=10, computeMarginalEffects=TRUE, showProgress=FALSE,
                                     method=vector("character", meth.length = ncol(x)), model=as.numeric(1)){

  y = yf10;  x = X;  method = "probit";  ndraw=1000 ;  burn.in=100;  thinning=1;  prior=list(a1=1, a2=1, beta=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 0);
  rho = 0.39;  rhosd = 0.16;  m=10;  computeMarginalEffects=TRUE;  showProgress=TRUE;
  
  if(!is.data.frame(y)){y = as.data.frame(y)}
  if(!is.matrix(x)){x = as.matrix(x)} 
  
  Yk = ncol(y)
  if(Yk != length(method)){
    stop("definition of method needs to occur for every covariate, so if Y1 is continuous and Y2 is dichotomous; method = c('continuous', 'probit')")
  }
  
  obs = list(); Nobs = list(); Nmiss = list(); vmean = list(); 
  in.ones = list(); nmk = list(); missing = list(); y_imp = list(); #y_imp_old = list();
  S = list(); H = list(); mu = list(); trW.i = list()
  betadraws = list(); ind = list(); ind2 = list(); zmean = list();
  beta = list(); for(i in 1:ncol(y)){beta[[i]] = rep(0, (ncol(X) + (i-1)))}
  total = list(); direct = list(); indirect = list(); lower = list(); upper = list()
  
  if("orderedprobit" %in% method){rho = start$rho}
  if("probit" %in% method){rho = start$rho}
  
  for(i in 1:Yk){
    obs[[i]] = which(!is.na(y[,i]))
    missing[[i]] = which(is.na(y[,i]))
    Nobs[[i]] = length(y[!is.na(y[,i]),i])
    Nmiss[[i]] = length(y[is.na(y[,i]),i])
    y_imp[[i]] = matrix(0, ndraw*thinning, Nmiss[[i]])
    #y_imp_old[[i]] = matrix(0, ndraw*thinning, Nmiss[[i]]) Y_imp_old is functional to compare the cut model against the former imputation method
    
    # partition data
    assign(paste0("W", i), W[obs[[i]], obs[[i]]])
    assign(paste0("y", i), y[obs[[i]],i])
    assign(paste0("Wy", i), get(paste0("W", i)) %*% get(paste0("y", i)))
  }
}
   


length(which(as.vector(W1) != 0))
