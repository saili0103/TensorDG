library(mvtnorm)
library(rTensor)
library(MASS)
library(CVXR)
#library(plotly)
#setwd("/Users/saili/Nutstore Files/Nutstore/Domain generalization-Multi-indices")

source('~/dg-main-functions.R')

#generate the observed groups
q=2
n=300
p.vec<-c(60,8,8)
Niter=500

###Exp1: fixing a, w, varying r####
A.list=list()
aa=6; ww=5
for(rr in 2:4){
  set.seed(aa*ww*rr)
  TLE.out<-NULL #mean & max &median errors of the observed grps, mean & max errors of the test grps, mean errors of all the grps
  Meta.out<-NULL
  Meta.out2<-NULL
  OLS.out<-NULL
    r.vec<- c(rr*2,rr,rr)
    tsr.gen<-tensor.gen(p.vec,r.vec)
    coef.tsr=tsr.gen$coef.tsr
    #svd(k_unfold(coef.tsr,1)@data)$d[1:(r.vec[1]+1)]
    #svd(k_unfold(coef.tsr,2)@data)$d[1:(r.vec[2]+1)]
    #svd(k_unfold(coef.tsr,3)@data)$d[1:(r.vec[3]+1)]
    A.list[[1]]=1:aa; A.list[[2]]=1:aa # arm set
    Omega=expand.grid(1:ww,1:ww) #body groups
    A1<-cbind(rep(A.list[[2]],each=p.vec[3]),rep(1:p.vec[3],length(A.list[[2]]))) 
    A2<-cbind(rep(1:p.vec[2],length(A.list[[1]])),rep(A.list[[1]],each=p.vec[2]))
    Ob=cbind(rep(A.list[[2]],each=p.vec[3]),rep(1:p.vec[3],length(A.list[[2]]))) #observed groups
    Ob=rbind(Ob, cbind(rep(1:p.vec[2],length(A.list[[1]])),rep(A.list[[1]],each=p.vec[2])))
    Ob<-unique(rbind(Ob, as.matrix(Omega)))
    Omega1<-unique(Omega[,1]); Omega2<-unique(Omega[,2])
    coef.Omega<-coef.tsr@data[,Omega1,Omega2]
    coef.Omega<-as.tensor(coef.Omega)
    grps<-expand.grid(1:p.vec[2],1:p.vec[3])
    grps.test<-grps[is.na(match(data.frame(t(grps)), data.frame(t(Ob)))),]
    #generate data
    for(it in 1:Niter){
      ####generate observed data####
      y=list(); X=list();y.til=list(); X.til=list()
      y.all=list();X.all=list()
      for(i in 1:nrow(Ob)){
        X[[i]]=rmvnorm(n/2,rep(0,p.vec[1]),diag(1,p.vec[1]))
        y[[i]]=X[[i]]%*%coef.tsr@data[,Ob[i,1],Ob[i,2]]+rnorm(n/2,0,1)
        X.til[[i]]=rmvnorm(n/2,rep(0,p.vec[1]),diag(1,p.vec[1]))
        y.til[[i]]=X.til[[i]]%*%coef.tsr@data[,Ob[i,1],Ob[i,2]]+rnorm(n/2,0,1)
        X.all[[i]]=rbind(X[[i]],X.til[[i]])
        y.all[[i]]=c(y[[i]],y.til[[i]])
      }
      #generate training data
      test.X=list(); test.y=list();sse.test3=rep(NA,nrow(grps.test))
      for(i in 1:nrow(grps.test)){
        test.X[[i]]=rmvnorm(n/2,rep(0,p.vec[1]),diag(1,p.vec[1]))
        test.y[[i]]=test.X[[i]]%*%coef.tsr@data[,grps.test[i,1],grps.test[i,2]]+rnorm(n/2,0,1)
        OLS.beta<-lm(test.y[[i]]~test.X[[i]]-1)$coef
        sse.test3[i]<-sum((OLS.beta-coef.tsr@data[,grps.test[i,1],grps.test[i,2]])^2)
      }
      OLS.out <- c(OLS.out, sse.test3)
      ####Method1:TensorTL w/ DS####
      #step 0 compute OLS for the observed groups
      #OLS.list=OLS.compute(X,y,Ob, p.vec)
      #OLS.list2=OLS.compute(X.til,y.til,Ob, p.vec)
      OLS.list.all=OLS.compute(X.all,y.all,Ob, p.vec)
      #step 1
      V0.hat<-SVD.0(OLS.list=OLS.list.all, Ob)
      SVD.re<-SVD.t(OLS.list=OLS.list.all, Ob, Omega=Omega, A.list=A.list)
      V.list=SVD.re$Vt

      #step 2
      B0.hat<-k_unfold(OLS.list.all$OLS.tsr,1)@data
      R0.hat<-V0.hat%*%solve(t(V0.hat)%*%B0.hat%*%t(B0.hat)%*%V0.hat)%*%t(V0.hat)%*%B0.hat%*%t(B0.hat)
      R1.hat<-Rt.compute(OLS.tsr=OLS.list.all$OLS.tsr, Ob=Ob,t=1, At=A.list[[1]],Omegat=Omega1, 
                         Vt.til=V.list[[1]], Til.tsr=OLS.list.all$OLS.tsr)
      R2.hat<-Rt.compute(OLS.tsr=OLS.list.all$OLS.tsr, Ob=Ob,t=2, At=A.list[[2]],Omegat=Omega2, 
                         Vt.til=V.list[[2]], Til.tsr=OLS.list.all$OLS.tsr)
      #step 3
      Beta.Omega.tsr<-as.tensor(OLS.list.all$OLS.tsr@data[,Omega1,Omega2])
      Beta.Omega0.tsr<-ttm(Beta.Omega.tsr,R0.hat,1) 
      Beta.hat.tsr<-ttm(Beta.Omega0.tsr,t(R1.hat),2)
      Beta.hat.tsr<-ttm(Beta.hat.tsr,t(R2.hat),3)
      
      sse.tl.out=rep(0,nrow(grps.test))
      for(j in 1:nrow(grps.test)){
        beta.tl=TL.dg2(X.tar=test.X[[j]],y.tar=test.y[[j]], grp.ind=as.numeric(grps.test[j,]), Beta.hat.tsr)
        sse.tl.out[j]=sum((beta.tl-coef.tsr@data[,grps.test[j,1],grps.test[j,2]])^2)
      }
      TLE.out<-c(TLE.out, sse.tl.out)
      ####Method 2: Meta-LM ####
      M.hat=diag(0,p.vec[1])
      for(i in 1:nrow(Ob)){
        for(k in 1:n){
          M.hat=M.hat+y.all[[i]][k]^2*X.all[[i]][k,]%*%t(X.all[[i]][k,])/n
        }
      }
      M.hat=M.hat/nrow(Ob)
      B.hat<-svd(M.hat)$u[,1:r.vec[1]]
      sse.meta.out=rep(0,nrow(grps.test))
      for(j in 1:nrow(grps.test)){
        alpha.hat<-lm(test.y[[j]]~test.X[[j]]%*%B.hat-1)$coef
        beta.meta=B.hat%*%alpha.hat
        sse.meta.out[j]=sum((beta.meta-coef.tsr@data[,grps.test[j,1],grps.test[j,2]])^2)
      }
      Meta.out<-c(Meta.out, sse.meta.out)
      ####Method 3: Meta-LM* ####
      B.hat=V0.hat
      sse.meta.out=rep(0,nrow(grps.test))
      for(j in 1:nrow(grps.test)){
        alpha.hat<-lm(test.y[[j]]~test.X[[j]]%*%B.hat-1)$coef
        beta.meta=B.hat%*%alpha.hat
        sse.meta.out[j]=sum((beta.meta-coef.tsr@data[,grps.test[j,1],grps.test[j,2]])^2)
      }
      Meta.out2<-c(Meta.out2, sse.meta.out)

     }
        
    write.table(TLE.out, file=paste('TL-TL-Exp1-s0',rr,'-',aa,'-',ww,'.txt',sep=""), col.names=F,row.names=F)
    write.table(Meta.out, file=paste('TL-Meta-Exp1-s0', rr,'-',aa,'-',ww,'.txt',sep=""), col.names=F,row.names=F)
    write.table(Meta.out2, file=paste('TL-Meta2-Exp1-s0', rr,'-',aa,'-',ww,'.txt',sep=""), col.names=F,row.names=F)
    write.table(OLS.out, file=paste('TL-OLS-Exp1-s0', rr,'-',aa,'-',ww,'.txt',sep=""), col.names=F,row.names=F)
    

  }#rr






