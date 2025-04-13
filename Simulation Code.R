library(mvtnorm)
library(MASS)
library(lme4)
library(ncvreg)
library(survival)
library(abind)
library(pracma)
library(pcoxtime)
library(dplyr)
library(lmerTest)
library(nlme)##lme
#options(warn=-1)
T=1
nrep=10
n=300  #sample size
p=5000# the dimension of mediators
S=1:9
q=2  #the dimenson of covariates
alpha=c(c(0.5,0.4,0.3),c(0.5,0.4,0.3),rep(0,3),rep(0,p-9))
eta=matrix(0.4,q,p)
theta1=0.4
theta2=rep(0.4,q)
beta=c(c(0.5,0.4,0.3),rep(0,3),c(0.5,0.4,0.3),rep(0,p-9))
true=1:3
S1=which(beta!=0)
for(irep in 1:nrep){
  ii=(T-1)*nrep+irep
  set.seed((T-1)*nrep+irep+1234)
  a1=rnorm(n,0,0.2)
  J1=matrix(rep(1:q,times=q),byrow=FALSE, nrow=q)
  K1=matrix(rep(1:q,times=q),byrow=TRUE, nrow=q)
  sigma=0.2^abs(J1-K1)
  Z=mvrnorm(n,rep(0,q),sigma)
  X=rnorm(n,0,0.5)
  # Consider intervals D1 = [0,t1), D2 = [t1,t2), D3 = [t2,t3), D4 = [t3,+infty)
  interval.max = 2
  delta1=0.2##20% center rate ,delta=0.2, 60% censer rate
  t=0:(interval.max-1)*delta1
  # Coefficient vector
  beta0.ti=c(theta1,theta2)
  beta0.tv=beta
  beta0 = c(beta0.ti, rep(beta0.tv, each = interval.max))
  p.ti=1+q
  p.tv=p
  p.x = p.ti + p.tv*interval.max
  M=matrix(0,n,p*interval.max)
  for(k in 1:p){
    b1=rnorm(n,0,0.2)
    for(j in 1:interval.max){
      M[,(k-1)*interval.max+j]=as.vector(X*alpha[k]+Z%*%eta[,k])+a1+b1+rnorm(n,0,1)
    }
  }
  ########## Austin Weibull #################
  delta2=1
  scale=1
  nu=4
  xx0 = cbind(X,Z,M,1:n)
  u = runif(n, 0, 1)
  # Calculate R1-R4, H1-H4
  H = R = Hinv = flag = matrix(0, nrow = n, ncol = interval.max)
  for (i in 1:interval.max){
    ind = c(1:p.ti, p.ti+(0:(p.tv-1))*interval.max+i)
    H[,i] = scale*exp(xx0[,ind]%*%c(beta0.ti,beta0.tv))
    if (i>1){
      R[,i] = R[,(i-1)] + H[,(i-1)]*(t[i]^nu-t[i-1]^nu)
      flag[,i] = -log(u)<R[,i]
    }
    Hinv[,i] = ((-log(u)-R[,i])/H[,i]+t[i]^nu)^(1/nu)
    Hinv[,i][which(Hinv[,i]=="NaN")]=Inf
  }
  h.ind = interval.max-rowSums(flag)
  for (i in 1:interval.max){
    ind.col = c(1:p.ti, p.ti+(0:(p.tv-1))*interval.max+i, dim(xx0)[2])
    ind.elig = h.ind>=i
    t0 = t[i]
    t1 = Hinv[ind.elig,i]
    if (i==1){
      t1[t1>t[i]+delta1] = t[i]+delta1
      status = t1==Hinv[ind.elig,i]
      xx = cbind(t0,t1,status,xx0[ind.elig, ind.col])
    }else{
      t1[t1>t[i]+delta2] = t[i]+delta2
      status = t1==Hinv[ind.elig,i]
      xx = rbind(xx, cbind(t0,t1,status,xx0[ind.elig, ind.col]))
    }
  }
  dim(xx)
  data=xx[order(xx[,ncol(xx)]),]
  rate=1-length(which(data[,"status"]!=0))/n
  rate
  id1=data[,ncol(data)]
  aa=table(id1)
  id2=unlist(lapply(1:length(aa), function(x)c(1:as.numeric(aa)[x])))
  tstart=data[,"t0"]
  tstop=data[,"t1"]
  status=data[,"status"]
  Xt=data[,4]
  Zt=data[,5:(4+q)]
  colnames(Zt)=1:q
  Mt=data[,(4+q+1):(4+q+p)]
  simData=data.frame(data,id2=id2)
  ############### SIS  step########################
  abhat=rep(0,p)
  count=rep(0,p)
  count_SIS=rep(0,p)
  count_naive=rep(0,p)
  beta.hat=rep(0,p)
  alpha.hat=rep(0,p)
  IDhat=integer(0)
  beta_SIS=beta_all=rep(0,p)
  alpha_SIS=alpha_pvalue=rep(0,p)
  for(k in 1:p){
    #print(k)
    fit1=coxph(Surv(tstart,tstop,status)~Xt+Zt+Mt[,k],cluster=id1)
    beta_SIS[k]=coef(fit1)[q+2]
    beta_all[k]=summary(fit1)$coef[q+2,6]
    Mk=Mt[,k]
    tryCatch({
      fit2=lme(Mk~Xt+Zt,random=~1|id1,control = lmeControl(opt = "optim"))
      alpha_SIS[k]=summary(fit2)$tTable[2,1]
      alpha_pvalue[k]=summary(fit2)$tTable[2,5]
    }, error = function(e) {
      cat("error", "\n")
      alpha_SIS[k]<<-0
      alpha_pvalue[k]<<-1
    })
  }
  d_0=round(n/(2*log(n)))
  ab_SIS=alpha_SIS*beta_SIS
  ID_SIS=which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d_0])
  count_SIS[ID_SIS]=1
  ############lasso for selecting beta#######################
  cv_fit=pcoxtimecv(Surv(tstart,tstop,status)~Xt+Zt+Mt[,ID_SIS],data=simData)
  print(cv_fit$lambda.min)
  fit3=pcoxtime(Surv(tstart,tstop,status)~Xt+Zt+Mt[,ID_SIS],data=simData,lambda=cv_fit$lambda.min)
  ID_lasso=ID_SIS[which(coef(fit3)[-c(1:(q+1))]!=0)]
  fit4=coxph(Surv(tstart,tstop,status)~Xt+Zt+Mt[,ID_lasso])
  betahat=summary(fit4)$coef[-c(1:(q+1)),1]
  beta_pvalue=summary(fit4)$coef[-c(1:(q+1)),5]
  d=length(ID_lasso)
  ################the multiple-testing  procedure ####
  bind.p <-rbind(pmin(beta_pvalue*d,1), pmin(alpha_pvalue[ID_lasso]*d,1))
  final.p <- apply(bind.p, 2, max)
  ID_fdr=ID_lasso[which(final.p<=0.05)]
  if(length(ID_fdr) > 0){
    fit6=coxph(Surv(tstart,tstop,status)~Xt+Zt+Mt[,ID_fdr])
    betahat=coef(fit6)[-c(1:(q+1))]
    beta.hat[ID_fdr]=betahat
    IDhat=ID_fdr
    alpha.hat=alpha_SIS
    abhat=beta.hat*alpha.hat
    count[ID_fdr]=1
  }
  save_res=list(IDhat=IDhat,abhat=abhat,beta.hat=beta.hat,alpha.hat=alpha.hat,count=count,count_SIS=count_SIS,
                count_naive=count_naive,abhat_naive=abhat_naive)
  save_name = paste('~/intermediate_results',ii,'.Rdata',sep="")
  save(save_res,file = save_name)
}
