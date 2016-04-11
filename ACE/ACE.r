library(mvtnorm)
library(splines)
library(truncnorm)

n=500; n1=250; p=2


#Total and burn-in iterations
MCAX=200; GNUM=180

Y=matrix(0,nrow=n,ncol=p)
X=array(0,c(n,2,p))
A=matrix(c(1,0.5,0.5,1),nrow=2); InvA=solve(A);

Ua=matrix(0,nrow=n,ncol=2)
Uc=matrix(0,nrow=n,ncol=2)

# Beta=matrix(0,nrow=MCAX,ncol=2)
Beta=c(1,0.5)
Sigma_A=0.6
Sigma_C=0.2
Sigma_E=0.1

for (i in 1:n1) {
    Ua[i,1]=rnorm(1,0,sqrt(Sigma_A)); Ua[i,2]=Ua[i,1]
    Uc[i,1]=rnorm(1,0,sqrt(Sigma_C)); Uc[i,2]=Uc[i,1]
    Ue=rnorm(2,0,sqrt(Sigma_E))
    #X[i,1:2,1]=rep(runif(1,5,13),2)
    X[i,1:2,1]=rnorm(2,0,2)
    X[i,1:2,2]=rnorm(2,0,2)
    Y[i,]=X[i,,]%*%Beta+Ua[i,]+Uc[i,]+Ue
}

for (i in (n1+1):n) {
    Ua[i,]=rmvnorm(1,c(0,0),Sigma_A*A)
    Uc[i,1]=rnorm(1,0,sqrt(Sigma_C)); Uc[i,2]=Uc[i,1]
    Ue=rnorm(2,0,sqrt(Sigma_E))
    #X[i,1:2,1]=rep(runif(1,5,13),2)
    X[i,1:2,1]=rnorm(2,0,2)
    X[i,1:2,2]=rnorm(2,0,2)
    Y[i,]=X[i,,]%*%Beta+Ua[i,]+Uc[i,]+Ue   ###This is y*
}

ua<-Ua
uc<-Uc

Ystar=Y
Y=exp(Ystar)

##Genrate Y##

kk=25

##Calculation of Penalty Matrix: M
D0=matrix(0,nrow=kk-1,ncol=kk)
D1=matrix(0,nrow=kk-1-1,ncol=kk-1)
d=c(-1,1)
for(dd in 1:(kk-1)){ D0[dd,dd:(dd+1)]=d }
for(dd in 1:(kk-2)){ D1[dd,dd:(dd+1)]=d }
M=t(D1%*%D0)%*%(D1%*%D0)
####

##B-splines basis functions and their first derivatives
Bk=array(0,c(kk,n,2))
bk=array(0,c(kk,n,2))

ks=matrix(nrow=2,ncol=kk+4) #knots

for (j in 1:2) {
    ks[j,1:4]=rep(min(Y[,j]),4)
    ks[j,(kk+1):(kk+4)]=rep(max(Y[,j]),4)
    
    for (k in 5:kk) { ks[j,k]=min(Y[,j])+(k-4)*(max(Y[,j])-min(Y[,j]))/(kk-3) }
    
    Bk[,,j]=t(spline.des(knots=ks[j,],x=Y[,j],ord=4)$design)
    bk[,,j]=t(spline.des(knots=ks[j,],x=Y[,j],ord=4,derivs=rep(1,n))$design)
}
###

#Initial values: Gam and tao
# ga1=seq(from=-1,to=8,length=kk)
# ga2=seq(from=-1,to=7,length=kk)
Gam=matrix(0,ncol=2,nrow=kk)
Gam[,1]<-seq(from=0.285,to=13.5,length=kk)
Gam[,2]<-seq(from=0.285,to=15,length=kk)
tao=1/rgamma(2,shape=1,rate=0.005)

###

#Tuning parameters and acceptance rate
sigma_r=c(2.7,2.7)
p_accept=rep(0,2)   
###

#Aid matrix
f=matrix(0,nrow=2,ncol=10)
###

#Grid points and their B-splines values
y<-matrix(0,nrow=2,ncol=100)
for (i in 1:2) {y[i,]=seq(from=min(Y[,i]),to=max(Y[,i]),length.out=100)}
Bk1=t(spline.des(knots=ks[1,],x=y[1,],ord=4)$design)
Bk2=t(spline.des(knots=ks[2,],x=y[2,],ord=4)$design)
BK=array(0,c(kk,100,2))
BK[,,1]<-Bk1
BK[,,2]<-Bk2

#Estimates of g(y)
GY<-array(0,c(100,MCAX-GNUM,2))


# Beta
# Beta=c(0,0)


gy<-matrix(0,nrow=n,ncol=p)

MCMC=matrix(0,nrow=MCAX,ncol=5)

Sigma_A=0.5
Sigma_C=0.5
for (i in 1:n1) {
    Ua[i,1]=rnorm(1,0,sqrt(Sigma_A)); Ua[i,2]=Ua[i,1]
    Uc[i,1]=rnorm(1,0,sqrt(Sigma_C)); Uc[i,2]=Uc[i,1]
}
 
for (i in (n1+1):n) {
    Ua[i,]=rmvnorm(1,c(0,0),Sigma_A*A)
    Uc[i,1]=rnorm(1,0,sqrt(Sigma_C)); Uc[i,2]=Uc[i,1]
}


#compare ystar and gy
com_gy=array(0,c(n,2,10))
com_gy_beta=array(0,c(n,2,10))
y_beta=matrix(0,nrow=n,ncol=p)


###gibs sampling
for (GIB in 1:MCAX) {
    
	for (j in 1:2) {
	
		Invsigma_e=Bk[,,j]%*%t(Bk[,,j])/Sigma_E+M/tao[j]   
		sigma_e=solve(Invsigma_e)
    
		ek=rmvnorm(1,sigma=sigma_e)
		unit_e=ek/sqrt(sum(ek*ek)) #standardize e
    
		lo_up=matrix(NA,ncol=kk,nrow=2)
        lo_up[1,]=rep(-Inf,kk)
		lo_up[2,]=rep(Inf,kk)
    
		for(ii in 1:(kk-1)){
			if (unit_e[ii+1]>unit_e[ii]) {
            lo_up[1,ii]=(Gam[ii,j]-Gam[ii+1,j])/(unit_e[ii+1]-unit_e[ii])
			}
			else {
            lo_up[2,ii]=(Gam[ii,j]-Gam[ii+1,j])/(unit_e[ii+1]-unit_e[ii])
			}
		}
    
		lower=max(lo_up[1,],na.rm=T)
		upper=min(lo_up[2,],na.rm=T)


    
		r1=rtruncnorm(10, a=lower, b=upper, mean=0, sd=sigma_r[j])
    
		w1=matrix(0,nrow=kk,ncol=10)
    
		for(ww in 1:10) { w1[,ww]=Gam[,j]+r1[ww]*unit_e }
    
		for(ff in 1:10) {
			f1=0
        
			for (i in 1:n) { 
            f1=f1-(t(w1[,ff])%*%Bk[,i,j]-Ua[i,j]-Uc[i,j])^2/(2*Sigma_E)+log(t(w1[,ff])%*%bk[,i,j])
			}
        
			f[1,ff]=f1-t(w1[,ff])%*%M%*%w1[,ff]/(2*tao[j])-(r1[ff])^2/(2*(sigma_r[j]^2))
        
		}
    
    

		f1_m=max(f[1,])
		
		for (ff in 1:10) { f[1,ff]=exp(f[1,ff]-f1_m) }
    
		p1=f[1,]/sum(f[1,])
    
		m=sample(1:10,1,prob=p1)  ##Sample 1 from 1 to 10 with probabilities: p1
		wstar=w1[,m]
    
    #Adding part: (9.1)
		lo_up=matrix(NA,ncol=kk,nrow=2)
        lo_up[1,]=rep(-Inf,kk)
        lo_up[2,]=rep(Inf,kk)
    
		for(ii in 1:(kk-1)){
			if (unit_e[ii+1]>unit_e[ii]) {
            lo_up[1,ii]=(wstar[ii]-wstar[ii+1])/(unit_e[ii+1]-unit_e[ii])
			}
			else {
            lo_up[2,ii]=(wstar[ii]-wstar[ii+1])/(unit_e[ii+1]-unit_e[ii])
			}
		}
    
		lower=max(lo_up[1,],na.rm=T)
		upper=min(lo_up[2,],na.rm=T)

    ###    print(c(lower, upper))
    ###
    
		w2=matrix(0,nrow=kk,ncol=10)
		r2=rtruncnorm(9, a=lower, b=upper, mean=0, sd=sigma_r[j])
    
		for(ww in 1:9) { w2[,ww]=wstar+r2[ww]*unit_e }
    
		w2[,10]=Gam[,j]
		r2[10]=r1[m]
    
		for(ff in 1:10) {
			f2=0
        
			for (i in 1:n) {
            f2=f2-(t(w2[,ff])%*%Bk[,i,j]-Ua[i,j]-Uc[i,j])^2/(2*Sigma_E)+log(t(w2[,ff])%*%bk[,i,j])
			}
        
			f[2,ff]=exp(f2-t(w2[,ff])%*%M%*%w2[,ff]/(2*tao[j])-(r2[ff])^2/(2*(sigma_r[j]^2))-f1_m)
		}

 #       print(f)

		Accept=sum(f[1,])/sum(f[2,])
    
		u=runif(1)
    
		if (u<=Accept) { #Accept: Update Gam, add 1 to acceptance rate
			Gam[,j]=wstar
            print(Accept)
            print(wstar)
			p_accept[j]=p_accept[j]+1
		}
 #       print(sum(wstar<0))
    
		tao[j]=1/rgamma(1, shape=1+(kk-2)/2, rate=0.005+t(Gam[,j])%*%M%*%Gam[,j]/2)

        #update tao
    ## print(t(Gam[,j])%*%M%*%Gam[,j])
		#burn-in
		if (GIB>GNUM) {
			GY[,GIB-GNUM,j]=t(Gam[,j])%*%BK[,,j]
		}
	
    }
	#update gy
    for(i in 1:2){
        gy[,i]<-t(Gam[,i])%*%Bk[,,i]
    }
    # print(gy)
    
    #compare
    if (GIB>(MCAX-10)){
        for(i in 1:n){
            y_beta[i,]<-X[i,,]%*%Beta
        }
        for(i in 1:2) {
            com_gy[,i,(GIB-(MCAX-10))]<-t(Gam[,i])%*%Bk[,,i]
            com_gy_beta[,i,(GIB-(MCAX-10))]<-t(Gam[,i])%*%Bk[,,i]-y_beta[,i]
        }
    }

#update Ua
    Sigma_A2=1/(2/Sigma_E+1/Sigma_A)
    for (i in 1:n1) {
        Mu_A=sum(gy[i,]-X[i,,]%*%Beta-Uc[i,])*Sigma_A2/Sigma_E
        Ua[i,1]=rnorm(1,Mu_A,sqrt(Sigma_A2)); Ua[i,2]=Ua[i,1]
    }

    Sigma_A2=solve(diag(1/Sigma_E,2)+InvA/Sigma_A)
    for (i in (n1+1):n) {
        Mu_A=Sigma_A2%*%(gy[i,]-X[i,,]%*%Beta-Uc[i,])/Sigma_E
        Ua[i,]=rmvnorm(1,Mu_A,Sigma_A2)
    }

#update Uc
    Sigma_C2=1/(2/Sigma_E+1/Sigma_C)
    for (i in 1:n) {
        Mu_C=sum(gy[i,]-X[i,,]%*%Beta-Ua[i,])*Sigma_C2/Sigma_E
        Uc[i,1]=rnorm(1,Mu_C,sqrt(Sigma_C2)); Uc[i,2]=Uc[i,1]
    }

#update Sigma_A
    a=(n1+2*(n-n1)+18)/2; b=4
    b=b+sum(Ua[1:n1,1]^2)
      
    for (i in (n1+1):n) {
        b=b+as.double(t(Ua[i,])%*%InvA%*%Ua[i,])
    }
     
    b=b/2
    Sigma_A=1/rgamma(1,shape=a,rate=b);

    MCMC[GIB,1]=Sigma_A

# update Sigma_C

    a=(n+18)/2; b=sum(Uc[,1]^2)/2+4
    Sigma_C=1/rgamma(1,shape=a,rate=b)

    MCMC[GIB,2]=Sigma_C

# update Sigma_E

    a=(2*n+18)/2; b=0
    for (i in 1:n) {     #####改到这里
    temp=gy[i,]-X[i,,]%*%Beta-Ua[i,]-Uc[i,]
    b=b+sum(temp*temp)
    }

    b=b/2 +4
    Sigma_E=1/rgamma(1,shape=a,rate=b);

    MCMC[GIB,3]=Sigma_E

# update beta
    Sigma_beta=matrix(0,nrow=p,ncol=p)
    Mu_beta=matrix(0,nrow=p,ncol=1)
    for (i in 1:n) {
        Sigma_beta=Sigma_beta+t(X[i,,])%*%X[i,,]
        Mu_beta=Mu_beta+t(X[i,,])%*%(gy[i,]-Ua[i,]-Uc[i,])
    }
    Sigma_beta=solve(Sigma_beta)*Sigma_E
    Mu_beta=Sigma_beta%*%Mu_beta/Sigma_E
    Beta=t(rmvnorm(1,Mu_beta,Sigma_beta))
    MCMC[GIB,4:5]=Beta
    
}	


# #Acceptance rate	
# p_accept=p_accept/MCAX
# print(p_accept)

# #Plots
# gy_1=rowMeans(GY[,,1])
# logy=log(y[1,])

# plot(logy~y[1,],type="l",col="red",ylim=c(-5,8))
# points(gy_1~y[1,],type="l",lty=2)

# #Plots
# gy_2=rowMeans(GY[,,2])
# logy2=log(y[2,])

# plot(logy2~y[2,],type="l",col="red",ylim=c(-5,8))
# points(gy_2~y[2,],type="l",lty=2)



# #beta
# BETA<-colMeans(MCMC[(GNUM+1):MCAX,])[3:4]
# print(BETA)

# #sigma_A
# sigma_a<-mean(MCMC[(GNUM+1):MCAX,1])
# print(sigma_a)

# #sigma_C
# sigma_c<-mean(MCMC[(GNUM+1):MCAX,2])
# print(sigma_c)

# ##compare plot
# Ystar1<-matrix(0,nrow=n,ncol=p)
# for(i in 1:2){Ystar1[,i]<-sort(Ystar[,i])}
# Ystar_beta<-Ystar-y_beta
# Ystar_beta1<-matrix(0,nrow=n,ncol=p)
# for(i in 1:2){Ystar_beta1[,i]<-sort(Ystar_beta[,i])}

# com_gy1=array(0,c(n,2,10))
# com_gy_beta1=array(0,c(n,2,10))
# for(j in 1:10){
#     for(i in 1:2){
#         com_gy1[,i,j]<-sort(com_gy[,i,j])
#         com_gy_beta1[,i,j]<-sort(com_gy_beta[,i,j])
#     }
# }

# title_wd=c('1.png','2.png','3.png','4.png','5.png','6.png','7.png','8.png','9.png','10.png')
# title_name=c('1','2','3','4','5','6','7','8','9','10')
# ##plot gy-Ystar for 1
# for(i in 1:10){
# #     png(file=title_wd[i], bg="transparent")
#     plot(com_gy1[,1,i]~Ystar1[,1],type="l",col="red")
#     abline(a = 0,b = 1)
#     title(title_name[i])
# #     dev.off()
# }

# ##plot gy-Ystar for 2
# for(i in 1:10){
# #     png(file=title_wd[i], bg="transparent")
#     plot(com_gy1[,2,i]~Ystar1[,2],type="l",col="red")
#     abline(a = 0,b = 1)
#     title(title_name[i])
# #     dev.off()
# }

# ##plot gy-Ystar for 1 without beta
# for(i in 1:10){
# #     png(file=title_wd[i], bg="transparent")
#     plot(com_gy_beta1[,1,i]~Ystar_beta1[,1],type="l",col="red")
#     abline(a = 0,b = 1)
#     title(title_name[i])
# #     dev.off()
# }

# ##plot gy-Ystar for 2 without beta
# for(i in 1:10){
# #     png(file=title_wd[i], bg="transparent")
#     plot(com_gy_beta1[,2,i]~Ystar_beta1[,2],type="l",col="red")
#     abline(a = 0,b = 1)
#     title(title_name[i])
# #     dev.off()
# }
