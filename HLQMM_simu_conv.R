rm(list=ls())
source(file="fonctions.R")

dyn.load('lqmm.so')
library('ald')
library('lamW')
library(sp)
library(expm)
library(SparseGrid)
library(splancs)
library(mvtnorm)
library(raster)
library(rgeos)
library(rqpd)
library(Matrix)
library(VGAM)

  #setting the hyperparameters
  MC=1000
  n=6
  JI=c(10,100,1000,10000)
  tau=0.5
  rho=0.02
  p=8
  phiC=5
  psiC=2
  betaC=c(1,2,-2,-1,3,5.5,0.5,0.1)
  sigmaC=c(1,2,1,3,1,0.5)
  
  
  #random effects generation
  h=sample(c(0.33,0.66,0.9),size=2*n,replace=TRUE)
  VAR=matrix(NA,ncol=n,nrow=n)
  for (i in 1:n){
    for (j in 1:n){
      VAR[i,j]=1-h[i+j]^rho
      
    }
    VAR[i,i]=1
  }
  VAR=nearPD(VAR)$mat
  nuC=rnorm(n,mean=0,sd=sqrt(phiC))
  nuC=sqrtm(VAR)%*%nuC
  
  
  names=data.frame(t(c("variable","J","MC")))
  write.table(names,'resconvnormale.csv',append=TRUE,col.names=FALSE,row.names=FALSE)
  for (J in JI){
    print(J)
    
      
     
      for (m in 1:MC){
        
        #Data generation : can be switch to any generation respecting the assumptions from the paper 
        X=X=matrix(rnorm(p*n*J,mean=2,sd=2),ncol=p)
        #X=matrix(rlaplace(8*n*J,location=2,scale=3),ncol=8)
        # X=matrix(rbinom(8*n*J,size=1,p=0.33),ncol=8)
        
        #Response generation from the data
        epsilonC=c()
        epsilonN=c()
        for (i in 1:n){
          epsilonC=append(epsilonC,rALD(J,mu=0,sigma=sigmaC[i],p=tau))
        }
        
        YC=matrix(nrow=(n*J),ncol=1)
        Z=matrix(0,nrow=(n*J),ncol=n)
        for (i in 1:n){
          Z[((i-1)*J+1):(i*J),i]=1
          YC[((i-1)*J+1):(i*J)]=rep(psiC,J)+X[((i-1)*J+1):(i*J),]%*%betaC+Z[((i-1)*J+1):(i*J),]%*%nuC+epsilonC[((i-1)*J+1):(i*J)]
          
        }
        
        groups=c()
        for (compt in 1:n){
          groups=append(groups,rep(compt,times=J))
        }
        group=as.factor(groups)
        
        #the method
        resC=tryCatch(lqmmtrue(fixed=YC ~ X,random=~1,group=groups,covariance=VAR,tau=tau,data=data.frame(YC,X,groups),control=list(method="df",startQR=FALSE,verbose=FALSE,UP_tol=1e-04,UP_max_iter=100,LP_max_iter=700,LP_tol_ll=1e-3,LP_tol_theta=1e-3)),error=function(e){return(list(theta_x=rep(1,ncol(X)+1),scale=rep(1,n),phi=1,opt=list(low_loop=6,upp_loop=1)))})
        
        print(m)
        #writing results on a csv file
        write.table(data.frame(t(c(c("psi"),c(J),m,resC$theta_x[1])))  ,'resconvnormale.csv',append=TRUE,col.names=FALSE,row.names=FALSE)
        write.table(data.frame(t(c(c("beta"),c(J),m,resC$theta_x[2:9]))),'resconvnormale.csv',append=TRUE,col.names=FALSE,row.names=FALSE)
        write.table(data.frame(t(c(c("sigma"),c(J),m,resC$scale))) ,'resconvnormale.csv',append=TRUE,col.names=FALSE,row.names=FALSE)
        write.table(data.frame(t(c(c("phi"),c(J),m,resC$phi))),'resconvnormale.csv',append=TRUE,col.names=FALSE,row.names=FALSE)
        write.table(data.frame(t(c(c("low_loop"),c(J),m,resC$opt$low_loop))) ,'resconvnormale.csv',append=TRUE,col.names=FALSE,row.names=FALSE)
        write.table(data.frame(t(c(c("upp_loop"),c(J),m,resC$opt$upp_loop))),'resconvnormale.csv',append=TRUE,col.names=FALSE,row.names=FALSE)
      }}

