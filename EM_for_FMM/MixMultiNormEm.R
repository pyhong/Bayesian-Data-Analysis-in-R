####预设参数
library(mvtnorm)
library(rootSolve)
simulate<-function(){
  sigma_1<-diag(c(1,2))
  sigma_2<-diag(c(1.5,3))
  mu_1<-c(1,2)
  mu_2<-c(3,4)
  pi_1<-0.3
  pi_2<-1-pi_1
  N<-500
  pi<-c(pi_1,pi_2)
  
  MixMultiNorm<-function(x)
  {
    sum(apply(sp,1,function(s){log(x[1]*dmvnorm(s,x[2:3],diag(x[4:5]))+x[6]*dmvnorm(s,x[7:8],diag(x[9:10])))}))
  }
  
  ####仿真开始
  ####产生样本
  
  sp_1<-rmvnorm(N*pi_1,mu_1,sigma_1)
  sp_2<-rmvnorm(N*pi_2,mu_2,sigma_2)
  random_index<-sample(c(1:N))
  sp<-rbind(sp_1,sp_2)[random_index,]
  ####初始化估计值
  e_sigma_1<-diag(runif(2))
  e_sigma_2<-e_sigma_1+diag(runif(2))
  
  e_mu_1<-runif(2)
  e_mu_2<-e_mu_1+runif(2)
  
  e_pi_1<-runif(1)
  e_pi_2<-1-e_pi_1
  e_pi<-c(e_pi_1,e_pi_2)
  e_zij<-matrix(0,ncol=2,nrow=N)
  dMultiNorm<-matrix(0,ncol=N,nrow=2)
  
  epsilon<-1
  while(epsilon > 1e-6){
    old_mu_1<-e_mu_1
    old_mu_2<-e_mu_2
    old_sigma_1<-e_sigma_1
    old_sigma_2<-e_sigma_2
    old_pi<-e_pi
    ####E-step
    dMultiNorm[1,]<-dmvnorm(sp,e_mu_1,e_sigma_1)
    dMultiNorm[2,]<-dmvnorm(sp,e_mu_2,e_sigma_2)
    z<-e_pi%*%dMultiNorm
    for(i in 1:N)
      for(j in 1:2)
      {
        e_zij[i,j]<-e_pi[j]*dMultiNorm[j,i]/z[i]
      }
    
    ####M-step
    e_pi<-colSums(e_zij)/N
    e_mu_1<-e_zij[,1]%*%sp/e_pi[1]/N
    e_mu_2<-e_zij[,2]%*%sp/e_pi[2]/N
    dif_1<-apply(sp,1,function(x){x-e_mu_1})
    dif_2<-apply(sp,1,function(x){x-e_mu_2})
    e_sigma_1<-dif_1%*%diag(e_zij[,1])%*%t(dif_1)/e_pi[1]/N
    e_sigma_2<-dif_2%*%diag(e_zij[,2])%*%t(dif_2)/e_pi[2]/N
    epsilon<-sum((old_mu_1-e_mu_1)^2)+sum((old_mu_2-e_mu_2)^2)+sum((e_pi-old_pi)^2)
    +sum((old_sigma_1-e_sigma_1)^2)
    +sum((old_sigma_2-e_sigma_2)^2)
  }
  para<-c(e_pi[1],e_mu_1,diag(e_sigma_1),e_pi[2],e_mu_2,diag(e_sigma_2))
  loglikeli<-MixMultiNorm(para)
  hes<-hessian(MixMultiNorm,para,pert = 1e-4)
  se<-sqrt(diag(solve(-hes)))
  result<-list(pi=e_pi,mu_1=e_mu_1,mu_2=e_mu_2,sigma_1=e_sigma_1,sigma_2=e_sigma_2,loglikeli<-loglikeli,se=se)
  result
}


sim_times<-500

mu_1<-matrix(NA,ncol=2,nrow=sim_times)
mu_2<-matrix(NA,ncol=2,nrow=sim_times)
pi<-matrix(NA,ncol=2,nrow=sim_times)
sigma_1<-matrix(NA,ncol=2,nrow=sim_times)
sigma_2<-matrix(NA,ncol=2,nrow=sim_times)
se<-matrix(NA,ncol=10,nrow=sim_times)

for(times in 1:sim_times)
{
  result<-simulate()
  mu_1[times,]<-result$mu_1
  mu_2[times,]<-result$mu_2
  pi[times,]<-result$pi
  sigma_1[times,]<-diag(result$sigma_1)
  sigma_2[times,]<-diag(result$sigma_2)
  se[times,]<-result$se
}

index<-pi[,2]<0.8

sim_mean<-list(
  mu_1_mean=colMeans(mu_1[index,],na.rm=T),
  mu_2_mean=colMeans(mu_2[index,],na.rm=T),
  pi_mean=colMeans(pi[index,],na.rm=T),
  sigma_1_mean=colMeans(sigma_1[index,],na.rm=T),
  sigma_2_mean=colMeans(sigma_2[index,],na.rm=T),
  se=colMeans(se[index,],na.rm=T)
)

sim_sd<-list(
  mu_1_se=apply(mu_1[index,],2,function(x){sd(x,na.rm=T)}),
  mu_2_se=apply(mu_2[index,],2,function(x){sd(x,na.rm=T)}),
  pi_se=apply(pi[index,],2,function(x){sd(x,na.rm=T)}),
  sigma_1_se=apply(sigma_1[index,],2,function(x){sd(x,na.rm=T)}),
  sigma_2_se=apply(sigma_2[index,],2,function(x){sd(x,na.rm=T)}),
  get_see=get_see<-function(){
  c(.pi_se[1],mu_1[1],mu_1[2],sigma_1[1],sigma_1[2],pi_se[2],mu_2[1],mu_2[2],sigma_2[1],sigma_2[2])
  }
  
)

