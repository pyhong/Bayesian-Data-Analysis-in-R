###LMM的EM算法
tr <- function(x) sum(diag(x))
em.mixed <- function(y, x, z, beta, var0, var1,maxiter=2000,tolerance = 1e-0010)
{
  library(MASS)
  time <-proc.time()
  n <- nrow(y)
  q1 <- nrow(z)
  conv <- 1
  L0 <- loglike(y, x, z, beta, var0, var1)
  i<-0
  cat(" Iter. sigma0 sigma1 Likelihood",fill=T)
  repeat {
    if(i>maxiter) {conv<-0
                   break}
    V <- c(var1) * z %*% t(z) + c(var0) * diag(n)
    Vinv <- solve(V)
    xb <- x %*% beta
    resid <- (y-xb)
    temp1 <- Vinv %*% resid
    s0 <- c(var0)^2 * t(temp1)%*%temp1 + c(var0) * n - c(var0)^2 * tr(Vinv)
    s1 <- c(var1)^2 * t(temp1)%*%z%*%t(z)%*%temp1+ c(var1)*q1 -
      c(var1)^2 *tr(t(z)%*%Vinv%*%z)
    w <- xb + c(var0) * temp1
    var0 <- s0/n
    var1 <- s1/q1
    beta <- ginv( t(x) %*% x) %*% t(x)%*% w
    L1 <- loglike(y, x, z, beta, var0, var1)
    if(L1 < L0) { print("log-likelihood must increase, llikel <llikeO, break.")
                  conv <- 0
                  break
    }
    i <- i + 1
    cat(" ", i," ",var0," ",var1," ",L1,fill=T)
    if(abs(L1 - L0) < tolerance) {break} #check for convergence
    L0 <- L1
  }
  list(beta=beta, var0=var0,var1=var1,Loglikelihood=L0)
}

####计算对数似然函数
loglike<- function(y, x, z, beta, var0, var1)
{
  library(MASS)
  n<- nrow(y)
  V <- c(var1) * z %*% t(z) + c(var0) * diag(n)
  Vinv <- ginv(V)
  xb <- x %*% beta
  resid <- (y-xb)
  temp1 <- Vinv %*% resid
  (-.5)*( log(det(V)) + t(resid) %*% temp1 )
}

tr <- function(x) sum(diag(x))

###########仿真开始
library(MASS)
i<-200
j<-3
beta<-c(0.2,0.5)
x1=rnorm(i*j)
x<-cbind(rep(1,i*j),x1)
sigma_u<-1
sigma_e<-sqrt(0.5)
u<-rep(rnorm(i),each=j)
e<-rnorm(i*j)
y<-x%*%beta+u+e
beta<-lm(y~x1)$coefficient
z<-diag(1,i*j)
lm.result<-em.mixed(y, x, z, beta, sigma_e^2, sigma_u^2,maxiter=2000,tolerance = 1e-0010)
