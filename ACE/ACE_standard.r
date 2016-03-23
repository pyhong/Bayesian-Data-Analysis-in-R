ACE_standard_R<-function(){
library(mvtnorm)
n=500; n1=250; p=2
Y=matrix(0,nrow=n,ncol=p)
X=array(0,c(n,2,p))
A=matrix(c(1,0.5,0.5,1),nrow=2); InvA=solve(A);

Ua=matrix(0,nrow=n,ncol=2)
Uc = matrix(0, nrow = n, ncol = 2)


Est=matrix(0,10,5)
SEst=matrix(0,10,5)

for (CIR in 1:10) {
	Beta=c(1,0.5)
	Sigma_A=0.4
	Sigma_C=0.3
	Sigma_E=0.2

	for (i in 1:n1) {
		Ua[i,1]=rnorm(1,0,sqrt(Sigma_A)); Ua[i,2]=Ua[i,1]
		Uc[i,1]=rnorm(1,0,sqrt(Sigma_C)); Uc[i,2]=Uc[i,1]
		Ue=rnorm(2,0,sqrt(Sigma_E))
		X[i,1:2,1]=rep(runif(1,5,13),2)
		X[i,1:2,2]=rnorm(2,0,2)
		Y[i,]=X[i,,]%*%Beta+Ua[i,]+Uc[i,]+Ue
	}

	for (i in (n1+1):n) {
		Ua[i,]=rmvnorm(1,c(0,0),Sigma_A*A)
		Uc[i,1]=rnorm(1,0,sqrt(Sigma_C)); Uc[i,2]=Uc[i,1]
		Ue=rnorm(2,0,sqrt(Sigma_E))
		X[i,1:2,1]=rep(runif(1,5,13),2)
		X[i,1:2,2]=rnorm(2,0,2)
		Y[i,]=X[i,,]%*%Beta+Ua[i,]+Uc[i,]+Ue
	}

	MCMC=matrix(0,500,5)

	Beta=c(0,0);
	Sigma_A=Sigma_C=Sigma_E=1

	for (i in 1:n1) {
		Ua[i,1]=rnorm(1,0,sqrt(Sigma_A)); Ua[i,2]=Ua[i,1]
		Uc[i,1]=rnorm(1,0,sqrt(Sigma_C)); Uc[i,2]=Uc[i,1]
	}

	for (i in (n1+1):n) {
		Ua[i,]=rmvnorm(1,c(0,0),Sigma_A*A)
		Uc[i,1]=rnorm(1,0,sqrt(Sigma_C)); Uc[i,2]=Uc[i,1]
	}

	for (GIB in 1:500) {
		Sigma_A2=1/(2/Sigma_E+1/Sigma_A)
		for (i in 1:n1) {
			Mu_A=sum(Y[i,]-X[i,,]%*%Beta-Uc[i,])*Sigma_A2/Sigma_E
			Ua[i,1]=rnorm(1,Mu_A,sqrt(Sigma_A2)); Ua[i,2]=Ua[i,1]
		}
		
		Sigma_A2=solve(diag(1/Sigma_E,2)+InvA/Sigma_A)
		for (i in (n1+1):n) {
			Mu_A=Sigma_A2%*%(Y[i,]-X[i,,]%*%Beta-Uc[i,])/Sigma_E
			Ua[i,]=rmvnorm(1,Mu_A,Sigma_A2)
		}
		
		Sigma_C2=1/(2/Sigma_E+1/Sigma_C)
		for (i in 1:n) {
			Mu_C=sum(Y[i,]-X[i,,]%*%Beta-Ua[i,])*Sigma_C2/Sigma_E
			Uc[i,1]=rnorm(1,Mu_C,sqrt(Sigma_C2)); Uc[i,2]=Uc[i,1]
		}
		
		a=(n1+2*(n-n1)+1)/2; b=0
		b=b+sum(Ua[1:n1,1]^2)
		
		for (i in (n1+1):n) {
			b=b+as.double(t(Ua[i,])%*%InvA%*%Ua[i,])
		}
		
		b=b/2
		Sigma_A=1/rgamma(1,shape=a,rate=b);
		MCMC[GIB,1]=Sigma_A

		a=(n+1)/2; b=sum(Uc[,1]^2)/2
		Sigma_C=1/rgamma(1,shape=a,rate=b)
		MCMC[GIB,2]=Sigma_C
		
		a=(2*n+1)/2; b=0
		for (i in 1:n) {     #####改到这里
			temp=Y[i,]-X[i,,]%*%Beta-Ua[i,]-Uc[i,]
			b=b+sum(temp*temp)
		}
		
		b=b/2
		Sigma_E=1/rgamma(1,shape=a,rate=b);
		MCMC[GIB,3]=Sigma_E
		
		Sigma_beta=matrix(0,nrow=p,ncol=p)
		Mu_beta=matrix(0,nrow=p,ncol=1)
		
		for (i in 1:n) {
			Sigma_beta=Sigma_beta+t(X[i,,])%*%X[i,,]
			Mu_beta=Mu_beta+t(X[i,,])%*%(Y[i,]-Ua[i,]-Uc[i,])
		}
		Sigma_beta=solve(Sigma_beta)*Sigma_E
		Mu_beta=Sigma_beta%*%Mu_beta/Sigma_E
		Beta=t(rmvnorm(1,Mu_beta,Sigma_beta))
		MCMC[GIB,4:5]=Beta
	}

	Est[CIR,]=apply(MCMC[301:500,],2,mean)
	SEst[CIR,]=apply(MCMC[301:500,],2,sd)
}
  return(list(Est=Est,SEst = SEst))
}

