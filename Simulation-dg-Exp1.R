library(mvtnorm)
library(rTensor)
library(MASS)
library(CVXR)
#library(plotly)


source('~/dg-main-functions.R')

#generate the observed groups
q=2
n=300
p.vec<-c(60,8,8)
Niter=500

###Exp1: fixing a, w, varying r####
A.list=list()
re.dg<-NULL
re.ols<-NULL#OLS
re.tl<-NULL #TL result
re.mm<-NULL #Maximin result
aa=6; ww=5
for(rr in 2:4){
  set.seed(aa*ww*rr)
  DG.mat<-matrix(NA, nrow=Niter, ncol=7)
  OLS.mat<-matrix(NA, nrow=Niter, ncol=7)
  TLE.mat<-matrix(NA, nrow=Niter, ncol=7) #mean & max &median errors of the observed grps, mean & max errors of the test grps, mean errors of all the grps
  MM.mat<-matrix(NA, nrow=Niter, ncol=7)      #generate 
    r.vec<- c(rr*2,rr,rr)
    tsr.gen<-tensor.gen(p.vec,r.vec)
    coef.tsr=tsr.gen$coef.tsr
   # svd(k_unfold(coef.tsr,1)@data)$d[1:(r.vec[1]+1)]
   # svd(k_unfold(coef.tsr,2)@data)$d[1:(r.vec[2]+1)]
   # svd(k_unfold(coef.tsr,3)@data)$d[1:(r.vec[3]+1)]
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
        test.X[[i]]=rmvnorm(n,rep(0,p.vec[1]),diag(1,p.vec[1]))
        test.y[[i]]=test.X[[i]]%*%coef.tsr@data[,grps.test[i,1],grps.test[i,2]]+rnorm(n,0,1)
        OLS.beta<-lm(test.y[[i]]~test.X[[i]]-1)$coef
        sse.test3[i]<-sum((OLS.beta-coef.tsr@data[,grps.test[i,1],grps.test[i,2]])^2)
      }
      
      ####Method1:TensorDG w/ DS####
      #step 0 compute OLS for the observed groups
      OLS.list=OLS.compute(X,y,Ob, p.vec)
      OLS.list2=OLS.compute(X.til,y.til,Ob, p.vec)
      OLS.list.all=OLS.compute(X.all,y.all,Ob, p.vec)
      #step 1
      V0.hat<-SVD.0(OLS.list=OLS.list.all, Ob)
      SVD.re<-SVD.t(OLS.list=OLS.list.all, Ob, Omega=Omega, A.list=A.list)
      V.list=SVD.re$Vt
      #step 2
      #R0.hat<-V0.hat%*%t(V0.hat)
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
      ###compute estimation errors
      sse.ob<-GSSE(Beta.hat.tsr,coef.tsr, Ob)
      sse.Omega<-GSSE(Beta.hat.tsr,coef.tsr, Omega)
      sse.test<-GSSE(Beta.hat.tsr,coef.tsr, grps.test)
      DG.mat[it,]<-c(mean(sse.ob), max(sse.ob), median(sse.ob), mean(sse.test),max(sse.test),
                      median(sse.test), sum(c(sse.ob,sse.test))/nrow(grps))
      #TDplot(B.hat=Beta.Omega0.tsr, B.ora=coef.Omega)
      #Beta.Omega.tsr<-as.tensor(Beta.hat.tsr@data[,Omega1,Omega2])
      #TDplot(B.hat=Beta.Omega.tsr, B.ora=coef.Omega)
      
      ####Method 2: OLS estimate####
      sse.ob3<-GSSE(OLS.list.all$OLS.tsr,coef.tsr, Ob)
      OLS.mat[it,]<-c(mean(sse.ob3), max(sse.ob3), median(sse.ob3), mean(sse.test3),max(sse.test3),
                       median(sse.test3), sum(c(sse.ob3,sse.test3))/nrow(grps))
      

      ####Method 4:Maximin ####
      beta.mat=k_unfold(OLS.list.all$OLS.tsr,1)@data
      tar.grps=which(apply(beta.mat,2,function(x) sum(x^2)==0))
      beta.mat=beta.mat[,-tar.grps]
      Sigma=diag(1,p.vec[1])
      mm.re<-Maxmin(beta.mat,Sigma)
      beta.mm<-Maxmin(beta.mat,Sigma)$beta.mm
      
      sse.mm.out=sapply(1:nrow(grps.test), function(j) 
        sum((beta.mm-coef.tsr@data[,grps.test[j,1],grps.test[j,2]])^2))
      sse.mm.ob=sapply(1:nrow(Ob), function(j) 
        sum((beta.mm-coef.tsr@data[,Ob[j,1],Ob[j,2]])^2))
      rw.mm.ob=sapply(1:nrow(Ob), function(j)  #reward fn
        mean(y[[j]])^2-mean((y[[j]]-X[[j]]%*%beta.mm)^2))
      rw.ols.ob=sapply(1:nrow(Ob), function(j)  #reward fn
        mean(y[[j]])^2-mean((y[[j]]-X[[j]]%*%beta.mat[,j])^2))
      rw.dg.ob=sapply(1:nrow(Ob), function(j)  #reward fn
        mean(y[[j]])^2-mean((y[[j]]-X[[j]]%*%Beta.hat.tsr@data[,Ob[j,1],Ob[j,2]])^2))
      
      MM.mat[it,]<-c(mean(sse.mm.ob), max(sse.mm.ob) , mean(sse.mm.out),max(sse.mm.out), 
                     min(rw.dg.ob), min(rw.ols.ob),sum(c(sse.mm.ob, sse.mm.out))/nrow(grps))
    }
    write.table(DG.mat, file=paste('DG-Exp1',rr,'-',aa,'-',ww,'.txt',sep=""), col.names=F,row.names=F)
    write.table(OLS.mat, file=paste('OLS-Exp1', rr,'-',aa,'-',ww,'.txt',sep=""), col.names=F,row.names=F)
    write.table(MM.mat, file=paste('MM-Exp1', rr,'-',aa,'-',ww,'.txt',sep=""), col.names=F,row.names=F)
    
    re.dg<-rbind(re.dg, c(rr,aa, ww, colMeans(DG.mat))) #with TensorDG with DS
    re.ols<-rbind(re.ols, c(rr,aa, ww, colMeans(OLS.mat))) ##OLS
    re.mm<-rbind(re.mm, c(rr,aa, ww, colMeans(MM.mat))) ##Maximin
  }#rr



write.table(re.dg, file='DG-Exp1.txt', col.names=F,row.names=F)
write.table(re.ols, file='OLS-Exp1.txt', col.names=F,row.names=F)
#write.table(re.tl, file='TL-Exp2.txt', col.names=F,row.names=F)
write.table(re.mm, file='MM-Exp1.txt', col.names=F,row.names=F)


