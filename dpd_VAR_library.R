

## This version assumes zero mean
## Code for real data analysis with VAR(p) models
## n observations of d dim'l vectors are formatted by a matrix with "d" rows and "n" columns
##
## Required package: doParallel, vars



library(doParallel)
library(vars)

obj.f0.chol <-function(X,p,theta, alpha){ # X: observation matirx
                                   # p: the order of VAR model
                                   # theta: c(vec(c), vec(A_1),...,vec(A_p), vech(Sig))    
                                   
      d<-nrow(X); n<-ncol(X)
      
      out = esttoVARchol(theta, d, p);
      get.resi = VAR.sigma(X, p, out$A);
      Sig = out$S;
      Inv.Sig = solve(Sig);
      Resi = cbind(matrix(0,d,p), get.resi$Resi); 
      
      tmp = 0;
      if(alpha == 0){
        for( i in 1:n){
          tmp<- tmp + (Resi[,i]%*%Inv.Sig%*%Resi[,i])
        }
        value = n*log(det(Sig))+tmp
        
      } else{
        for( i in 1:n){
          tmp<-tmp + exp(-0.5*alpha*t(Resi[,i])%*%Inv.Sig%*%Resi[,i])
        }
        value = det(Sig)^(-alpha/2)*(n/(1+alpha)-(1+1/alpha)*tmp)
      }
   return(value)
}



obj.t0.chol <-function(obs,p,est,alpha){
  
  d<-nrow(obs); n = ncol(obs)
  #    est = c(rep(0, d), est);
  #    const<-est[1:d]
  out = esttoVARchol(est, d, p);
  get.resi = VAR.sigma(obs, p, out$A);
  Sig = out$S;
  Inv.Sig = solve(out$S);
  Resi = cbind(matrix(0,d,p), get.resi$Resi); 
  
  tmp = NULL;
  for( i in 1:n){
    tmp<-c(tmp,Resi[,i]%*%Inv.Sig%*%Resi[,i])
  }
  
  if(alpha==0){ log(2*pi)+0.5*log(det(Sig))+0.5*tmp
  } else {
    tmp<-exp(-0.5*alpha*tmp)
    (4*pi^2*det(Sig))^(-alpha/2)*(1/(1+alpha)-(1+1/alpha)*tmp)
  }
}


# fuction that returns mdpd estimates; alpha=0 yields ML estimates
mdpde0.chol <-function(obs,p,alpha,init=NULL){  
# if you have a initial parmeter vector,
# use the one as initial value  for optimization 
  
  obj<-function(theta){ obj.f0.chol(obs,p,theta,alpha)}
  d<-nrow(obs)
  
  if(is.null(init)){
#    const0<-rep(0,d);  A0<-rep(0,p*d^2)
#    Sig0<-rep(0.1, d*(d+1)/2)
#    k<-1; for(i in 1:d){Sig0[k]<-1; k<-k+d-i+1}
#    init<-c(const0,A0,Sig0)
#    ## Initial estimator from LSE
    aa = VAR.lse(obs, p=p);
    rr.sig = t(chol(aa$Sigma));
    init = c(as.vector(aa$hatA), rr.sig[lower.tri(rr.sig, diag = TRUE)]);
  }
#  print(init);
  # ff = constrOptim.nl(init,obj,hin=hin.g0,control.outer=list(trace=F))$par
  ff = optim(init, obj)$par;
  
  out = esttoVARchol(ff, d, p);
  out$est = ff;
  out$init = init;
  
  return(out);
}


## Parameter changes for A or Sigma or Both
score.test.chol.AS <-function (obs, p, est, alpha, type, D, n.cl = 4) 
{
  n <- ncol(obs)
  d = nrow(obs)
  if (missing(D)) {
    D = max(30, 0.05 * n)
  }
  if (missing(type)) {
    type = "AS"
  }
  if (missing(est)) {
    est0 = mdpde0.chol(obs, p = p, alpha = alpha, init = NULL)
    est = est0$est
  }
  if (type == "A") {
    n.est <- p * d * d
  }
  if (type == "S") {
    n.est <- d * (d + 1)/2
  }
  if (type == "AS") {
    n.est <- length(est)
  }
  get.ASscore.chol = function(subobs, p, est, alpha, type) {
    dx = 0.001
    d = nrow(subobs)
    n <- ncol(subobs)
    L0 <- obj.t0.chol(subobs, p, est, alpha)
    if (type == "A") {
      n.est <- p * d * d
      d.theta <- abs(est * dx + 1e-04)
      dh <- diag(d.theta[1:n.est])
      dh = rbind(dh, matrix(0, d * (d + 1)/2, n.est))
      score <- NULL
      for (j in 1:n.est) {
        score <- rbind(score, (obj.t0.chol(subobs, p, 
                                           est + dh[, j], alpha) - L0)/(dh[j, j]))
      }
    }
    if (type == "S") {
      n.est <- d * (d + 1)/2
      d.theta <- abs(est * dx + 1e-04)
      dh <- diag(d.theta[-(1:(p * d^2))])
      dh = rbind(matrix(0, p * d^2, n.est), dh)
      score <- NULL
      for (j in 1:n.est) {
        score <- rbind(score, (obj.t0.chol(subobs, p, 
                                           est + dh[, j], alpha) - L0)/(dh[(p * d^2 + 
                                                                              j), j]))
      }
    }
    if (type == "AS") {
      n.est <- length(est)
      d.theta <- abs(est * dx + 1e-04)
      dh <- diag(d.theta[1:n.est])
      score <- NULL
      cl <- makeCluster(n.cl)
      registerDoParallel(cl)
      if (n.cl > 1) {
        score = foreach(j = 1:n.est, .combine = rbind, 
                        .export = c("obj.t0.chol", "esttoVARchol", 
                                    "VAR.sigma")) %dopar% {
                                      (obj.t0.chol(subobs, p, est + dh[, j], alpha) - 
                                         L0)/(dh[j, j])
                                    }
        stopCluster(cl)
      }
      else {
        for (j in 1:n.est) {
          score <- rbind(score, (obj.t0.chol(subobs, 
                                             p, est + dh[, j], alpha) - L0)/(dh[j, j]))
        }
      }
    }
    hat.K <- matrix(0, n.est, n.est)
    for (i in 1:n) {
      hat.K <- hat.K + score[, i] %*% t(score[, i])
    }
    if (rcond(hat.K) > 0.001) {
      hat.K.Inv <- solve(hat.K) * n
      sigTF = FALSE
      temp <- rep(0, n)
      score.k <- rep(0, n.est)
      for (k in 1:n) {
        score.k <- score.k + score[, k]
        temp[k] <- t(score.k) %*% hat.K.Inv %*% score.k
      }
    }
    else {
      sigTF = TRUE
      temp = rep(0, n)
    }
    if (sum(is.nan(temp)) > 0) {
      sigTF = TRUE
    }
    return(list(tstat = temp, score = score, K = hat.K, 
                singularTF = sigTF))
  }
  temp = get.ASscore.chol(obs, p, est, alpha, type)
  time.index = seq(from = D + 1, to = n - D, by = 1)
  if (!temp$singularTF) {
    tscore0 <- rep(0, n)
    tscore0.k <- rep(0, n.est)
    I.inv = (solve(temp$K) * n)
    for (k in 1:n) {
      tscore0.k <- tscore0.k + temp$score[, k]
      tscore0[k] <- t(tscore0.k) %*% I.inv %*% tscore0.k
    }
    khat0 = D + which.max(tscore0[time.index])
    hat.Ic <- matrix(0, n.est, n.est)
    est.1 <- mdpde0.chol(obs[, 1:khat0], p = p, alpha = alpha, 
                         init = NULL)
    tmp1 <- get.ASscore.chol(obs[, (1:khat0)], p, est.1$est, 
                             alpha, type)
    est.2 <- mdpde0.chol(obs[, -(1:khat0)], p = p, alpha = alpha, 
                         init = NULL)
    tmp2 <- get.ASscore.chol(obs[, -(1:khat0)], p, est.2$est, 
                             alpha, type)
    hat.Ic <- (tmp1$K + tmp2$K)
    if (rcond(hat.Ic) < 0.001) {
      hat.Ic.Inv = I.inv
    }
    else {
      hat.Ic.Inv = solve(hat.Ic) * n
    }
    tscore <- rep(0, n)
    tscore.k <- rep(0, n.est)
    for (k in 1:n) {
      tscore.k <- tscore.k + temp$score[, k]
      tscore[k] <- t(tscore.k) %*% hat.Ic.Inv %*% tscore.k
    }
  }
  else {
    tscore = tscore0 = rep(0, n)
    khat0 = D + which.max(tscore0[time.index])
  }
  out = list()
  out$est = est
  out$score = tscore/n
  khat = D + which.max(tscore[time.index])
  out$maxscore = tscore[khat]/n
  out$khat = khat
  out$score0 = tscore0/n
  out$maxscore0 = (tscore0[khat0])/n
  out$khat0 = khat0
  return(out)
}




esttoVAR = function(est0, d, p){
  a = est0[1:(p*d^2)]; b= est0[-(1:(p*d^2))];
  hatA=matrix(a, byrow=FALSE, nrow=d);
  hatS = matrix(0, d, d);
  hatS[lower.tri(hatS, diag = TRUE)] = b;
  S = hatS + t(hatS);
  diag(S) = diag(hatS) 
  return(list(A = hatA, S=S))
}

esttoVARchol = function(est0, d, p){
  a = est0[1:(p*d^2)]; b= est0[-(1:(p*d^2))];
  hatA=matrix(a, byrow=FALSE, nrow=d);
  hatS = matrix(0, d, d);
  hatS[lower.tri(hatS, diag = TRUE)] = b;
  ## Use cholesky decompositon in parametrizaiton.
  S = hatS%*%t(hatS);
  return(list(A = hatA, S=S))
}

VAR.sigma = function (y, p, A) 
{
  k = dim(y)[1]
  Tt = dim(y)[2]
  T1 = Tt - p
  m = sum(A != 0)
  eps = 0.1
  Resi = mat.or.vec(k, T1)
  X = mat.or.vec(k * p, T1)
  Y1 = mat.or.vec(k, T1)
  for (j in 1:T1) {
    id = seq(from = j + p - 1, to = j, by = -1)
    x = as.vector(y[, id])
    Resi[, j] = y[, p + j] - A %*% x
  }
  Sigma_z = Resi %*% t(Resi)/T1
  LL = T1 * log(det(Sigma_z))
  BIC = LL + ((log(T1))^2) * (m + k * (k + 1)/2)
  out = list()
  out$Resi = Resi
  out$Sigma_z = Sigma_z
  out$BIC = BIC
  out$LL = LL
  out$m = m
  out$T1 = T1
  return(out)
}

################################
# VAR(p) best model 
# BIC is calculated based on RSS
################################

VAR.best = function(y, maxorder){
  ## Input data is dim*N
  N = ncol(y);
  dim = nrow(y);
  BIC = numeric(maxorder);
  for(p in 1:maxorder){
    ar1 = VAR(t(y), p=p, type="none")
    BIC[p] =  as.numeric(VAR.sigma(y, p=p, Bcoef(ar1))$BIC);
  }
  p.opt = which.min(BIC);
  best = VAR(t(y), p=p.opt, type="none")
  out = VAR.sigma(y, p=p.opt, Bcoef(best));
  out$hatA = Bcoef(best);
  out$BIC = BIC;
  out$opt = p.opt
  return(out)
}


###########################################
# OLS estimation of VAR
###########################################

VAR.lse  = function(y, p){
  y = y - rowMeans(y);
  Tt = dim(y)[2];
  T1 = Tt-p;
  k = dim(y)[1];
  
  # create vector X1 and Y1
  X1 = matrix(0,k*p, T1); 
  Y1 = matrix(0,k, T1);
  for(j in 1:T1){
    # ar term
    id = seq(from= j+p-1, to = j, by=-1);
    x = as.vector(y[,id]);
    X1[,j] = x;
    Y1[,j] = y[,(j+p)];
  }
  
  hatA = Y1%*%t(X1)%*%solve(X1%*%t(X1));
  
  Resi = matrix(0,k, T1); 
  # residuals
  for(j in 1:T1){
    id = seq(from= j+p-1, to = j, by=-1);
    x = as.vector(y[,id]);
    Resi[,j] =  y[,p+j]  - hatA%*%x;
  }
  Sigma_z = Resi%*%t(Resi)/T1;
  return(list(hatA = hatA, Sigma=Sigma_z, p=p))
}



### P-value / critical value calculation
# We need the quantiles of sup of ||26-dimensional standard Brownian brdige||_2^2.

q.Brown.bridge2<-function(D,n,r, level=.95, n.cl=4){
  
  if(missing(r)){ r = 2000;}
  dist<-rep(0,r);
  set.seed(5678)
  
  if(n.cl >1){
    cl <- makeCluster(n.cl)  
    registerDoParallel(cl) 
  dist= foreach(k=1:r, .combine=c) %dopar% {
    e <-matrix(rnorm(n*D), ncol=n); ## Generate D*n matrix
    B = apply(e, 1, cumsum); ## This is column
    m = rowMeans(e);
    BB = B - (1:n)%*%t(m);
    BB = BB/sqrt(n);
    max(apply(BB^2,1,sum))
  }
  stopCluster(cl);
  } else {  
  for(k in 1:r){
    e <-matrix(rnorm(n*D), ncol=n); ## Generate D*n matrix
    B = apply(e, 1, cumsum); ## This is column
    m = rowMeans(e);
    BB = B - (1:n)%*%t(m);
    BB = BB/sqrt(n);
    dist[k]<-max(apply(BB^2,1,sum))
  }
  }

  return(quantile(dist, level))  
}


pval.Brown.bridge<-function(crit, D, n, r){
  if(missing(r)){ r = 2000;}
  dist<-rep(0,r)
  for(k in 1:r){
    BB<-matrix(0,D,n)
    for(i in 1:D){
      e<-rnorm(n)
      m<-mean(e)
      BB[i,]<-(cumsum(e)-(1:n)*m )/sqrt(n)
    }
    dist[k]<-max(apply(BB^2,2,sum))
  }
  pval = sum(dist > crit)/r;
  return(pval)  
}

#################################################
# Binary segmentation with Cholesky decomposition
#################################################


score.chol.binary =function (obs, p, alpha = 0, type, D, n.cl, level = 0.95, critical, maxK = 3, 
                             plotTF = TRUE) 
{
  n = length = ncol(obs)
  d = dim = nrow(obs)
  if (missing(D)) {
    D = max(50, 0.05 * n)
  }
  if (missing(type)) {
    type = "AS"
  }
  if (missing(n.cl)) {
    n.cl = n.cl = detectCores(logical = TRUE) - 1
  }
  if (type == "A") {
    n.est <- p * d * d
  }
  if (type == "S") {
    n.est <- d * (d + 1)/2
  }
  if (type == "AS") {
    n.est <- p * d * d + d * (d + 1)/2
  }
  if(missing(critical)){ critical = q.Brown.bridge2(n.est, n = 1000, r = 1000, level = level, 
                             n.cl = n.cl);}
  out1 = score.test.chol.AS(obs, p = p, alpha = alpha, type = type, 
                            D = D, n.cl = n.cl)
  tstat = out1$maxscore
  khat = out1$khat
  result = list(tstathist = tstat, Brhist = khat)
  br = c(0, n)
  nmax = 1
  if (tstat > critical) {
    br = c(0, khat, n)
    st = sum(tstat > critical)
    spindex = c(1, 1)
    while (st > 0 && nmax < maxK) {
      nmax = nmax + 1
      lbr = length(br)
      Ts = NULL
      sst = Br = NULL
      brindex = seq(from = 1, to = lbr - 1, by = 1)
      for (j in brindex) {
        if (spindex[j]) {
          id = seq(from = br[j] + 1, to = br[j + 1], 
                   by = 1)
          dat = obs[, id]
          if (ncol(dat) > 2 * D + 10) {
            out2 = score.test.chol.AS(dat, p = p, alpha = alpha, 
                                      type = type, D = D, n.cl = n.cl)
            tstat2 = out2$maxscore
            khat2 = out2$khat
            Br = c(Br, br[j] + khat2)
            Ts = c(Ts, tstat2)
            idf = 1 * (tstat2 > critical)
            sst = c(sst, rep(idf, idf + 1))
          }
          else {
            sst = c(sst, 0)
          }
        }
        else {
          sst = c(sst, 0)
        }
      }
      st = sum(sst)
      if (!is.null(Ts)) {
        newbr = (abs(Ts) > critical) * Br
        br = unique(c(br, newbr))
      }
      br = sort(br)
      spindex = sst
      result$tstathist = c(result$tstathist, Ts)
      result$Brhist = c(result$Brhist, Br)
    }
  }
  if (plotTF) {
    plot(seq(1, length), out1$score, type = "l", xlab = "", 
         ylab = "", main = "Score")
    abline(v = br, col = "blue")
  }
  result$khat = br
  result$iter = length(br) - 1
  result$critical = critical
  result$score = out1$score
  return(result)
}





