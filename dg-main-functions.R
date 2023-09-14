library(glmnet)
library(CVXR)
#generating the true tensor model
tensor.gen<-function(p.vec,r.vec){
  #r.vec: length q+1 vector, rank
  #p.vec: length q+1 vector, dimension
  q<-length(p.vec)-1
  core.tsr <- rand_tensor(modes = r.vec)#core tensor
  tsr1<-k_unfold(core.tsr,1)
  tsr.uf1<-apply(tsr1@data,2, function(v) v/sqrt(sum(v^2)))
  core.tsr<-k_fold(tsr.uf1,m=1, modes=core.tsr@modes)
  R<-list()
  coef.tsr<-core.tsr
  for(i in 1:(q+1)){
    R[[i]]<-matrix(rnorm(r.vec[i]*p.vec[i]),ncol=p.vec[i]) 
    R[[i]][1:r.vec[i],1:r.vec[i]]=diag(1,r.vec[i])
    R[[i]]<-apply(R[[i]],2, function(v) v/sqrt(sum(v^2)))
    coef.tsr<- ttm(tnsr = coef.tsr, mat = t(R[[i]]), m = i)
  }
  return( list(coef.tsr=coef.tsr, R.list=R, core.tsr=core.tsr))
}

GSSE<-function(B.hat,B.ora, grps){
  sapply(1:nrow(grps), function(k) 
    sum((B.hat@data[,grps[k,1],grps[k,2]]-B.ora@data[,grps[k,1],grps[k,2]])^2))
}

tensor.decomp.ora<-function(coef.tsr, A.list, Omega, grps.test){
  #r.vec: length q+1 vector, rank
  #p.vec: length q+1 vector, dimension
  Omega1<-unique(Omega[,1])
  Omega2<-unique(Omega[,2])
  beta.Omega<-coef.tsr@data[,Omega1,Omega2]
  B0<-k_unfold(coef.tsr,1)@data
  V0<-svd(B0%*%t(B0))$u[,1:(2*rr)]
  B1.ar<-k_unfold(as.tensor(coef.tsr@data[,,A.list[[1]]]),2)@data
  B1.jo=B1.ar[Omega1,]
  V1<-svd(B1.jo%*%t(B1.jo))$u[,1:rr]
  svd(B1.jo%*%t(B1.jo))$d
  B2.ar<-k_unfold(as.tensor(coef.tsr@data[,A.list[[2]],]),3)@data
  B2.jo=B2.ar[Omega2,]
  V2<-svd(B2.jo%*%t(B2.jo))$u[,1:rr]
  
  Gam1<-solve(t(t(B1.jo)%*%V1)%*%t(B1.jo)%*%V1)%*%t(t(B1.jo)%*%V1)%*%t(B1.ar)
  Gam2<-solve(t(t(B2.jo)%*%V2)%*%t(B2.jo)%*%V2)%*%t(t(B2.jo)%*%V2)%*%t(B2.ar)
  Beta.Omega.ora<-as.tensor(beta.Omega)
  Beta.Omega0.tsr<-ttm(Beta.Omega.ora,V0%*%t(V0),1) 
  Beta.ora.tsr<-ttm(Beta.Omega0.tsr,t(V1%*%Gam1),2)
  Beta.ora.tsr<-ttm(Beta.ora.tsr,t(V2%*%Gam2),3)
  test.err<-mean(GSSE(Beta.ora.tsr, coef.tsr, grps.test))
  R.list=list(); V.list=list()
  R.list[[1]]<-V0%*%t(V0); R.list[[2]]<-V1%*%Gam1; R.list[[3]]<-V2%*%Gam2
  V.list[[1]]<-V0; V.list[[2]]<-V1; V.list[[3]]<-V2
  return(list(Beta.Omega.ora=Beta.Omega.ora, R.list=R.list, V.list=V.list, test.err=test.err))
}
OLS.compute<-function(X,y, Ob, p.vec){
  sig.mat<-matrix(0,p.vec[2],p.vec[3])
  n.mat<-matrix(0,p.vec[2],p.vec[3])
  invTr.mat<-matrix(0,p.vec[2],p.vec[3])
  bc0<-diag(0,p.vec[1])#bias-corrected term for mode 0
  OLS.arr<-array(0,dim=p.vec)
  OLS.arr2<-array(0,dim=p.vec)
  for(i in 1:nrow(Ob)){
    beta.hat<-lm(y[[i]]~X[[i]]-1)$coef
    OLS.arr[,Ob[i,1],Ob[i,2]]<-beta.hat
    sig.mat[Ob[i,1],Ob[i,2]]<-sum((y[[i]]-X[[i]]%*%beta.hat)^2)/(length(y[[i]])-p.vec[1])
    n.mat[Ob[i,1],Ob[i,2]]<-length(y[[i]])
    Sig.inv<-solve(cov(X[[i]]))
    invTr.mat[Ob[i,1],Ob[i,2]]<-sum(diag(Sig.inv))
    bc0<-bc0+sig.mat[Ob[i,1],Ob[i,2]]*Sig.inv/length(y[[i]])
  }
  list(OLS.tsr=as.tensor(OLS.arr), sig.mat=sig.mat,n.mat=n.mat,
                invTr.mat=invTr.mat,bc0=bc0) #for estimating Vt
}

#currently work for q=2
Rt.compute<-function(OLS.tsr, Ob, t, Omegat, At, Vt.til, rt=rt, Til.tsr=NULL){
  #t=1; At=A1; Omegat=Omega1; Vt.til=V.list[[1]]
  p=OLS.tsr@modes[1]
  pt=OLS.tsr@modes[t+1]
  if(is.null(Til.tsr)){
    Til.tsr=OLS.tsr
  }
  if(t==1){
    Beta.Ar.tsr<-as.tensor(OLS.tsr@data[,,At])
    Beta.jo.til<-k_unfold(as.tensor(Til.tsr@data[,Omegat,At]),t+1)@data
  }else{ #t=2
    Beta.Ar.tsr<-as.tensor(OLS.tsr@data[,At,])
    Beta.jo.til<-k_unfold(as.tensor(Til.tsr@data[,At,Omegat]),t+1)@data
  }
  Beta.Ar.mat<-k_unfold(Beta.Ar.tsr,t+1)@data
  Beta.jo.mat<-Beta.Ar.mat[Omegat,]
  BV.jo<-t(Beta.jo.mat)
  BV.jo.til<-t(Beta.jo.til)
  Gamt<-solve(t(BV.jo.til%*%Vt.til)%*%BV.jo%*%Vt.til)%*%t(BV.jo.til%*%Vt.til)%*%t(Beta.Ar.mat)
  #colMeans(t(Beta.Ar.mat)-BV.jo%*%Vt.til%*%Gamt)
  Rt<-Vt.til%*%Gamt
  Rt

}


SVD.t<-function(OLS.list, Ob, Omega, A.list,r.vec=NULL){
  OLS.tsr=OLS.list$OLS.tsr
  p=OLS.tsr@modes[1]
  bc0=OLS.list$bc0
  q=ncol(Omega)
  Omega1=unique(Omega[,1]);Omega2=unique(Omega[,2])
  C1=unique(c(Omega2,A.list[[1]])); C2=unique(c(Omega1,A.list[[2]]))
  Vt=list()
  for(t in 1:q){
    if(t==1){
      Beta.Omega.tsr<-as.tensor(OLS.tsr@data[,Omega1,C1])
    }else{
      Beta.Omega.tsr<-as.tensor(OLS.tsr@data[,C2,Omega2])
    }

    B.mat<-k_unfold(Beta.Omega.tsr, t+1)@data
    Omegat=unique(Omega[,t])
    tr.Sig<-rep(0,length(Omegat))
    for(j in 1:length(Omegat)){
      if(t==1){
        tr.Sig[j]=sum(OLS.list$sig.mat[Omega1[j],Omega2]*
                        OLS.list$invTr.mat[Omega1[j],Omega2]/OLS.list$n.mat[Omega1[j],Omega2])
        N= sum(OLS.list$n.mat[Omega1[j],Omega2])
      }else{
        tr.Sig[j]=sum(OLS.list$sig.mat[Omega1,Omega2[j]]*
                        OLS.list$invTr.mat[Omega1,Omega2[j]]/OLS.list$n.mat[Omega1,Omega2[j]])
        N= sum(OLS.list$n.mat[Omega1,Omega2[j]])
      }
    }
    Thetat.hat<-(B.mat)%*%t(B.mat)-diag(tr.Sig)
    eigent<-eigen(Thetat.hat)
    n.avg<-OLS.list$n.mat[Ob[1,1],Ob[1,2]]
    if(is.null(r.vec)){
      thres<-sqrt(eigent$values[1]*(length(Omegat)+log(n.avg))*length(A.list[[t]])/N)
      rt=sum(eigent$values>=thres)
    }else{rt<-r.vec[t+1]}
    Vt[[t]]<-eigent$vectors[,1:rt] 
  }
  return(list(Vt=Vt))
}


SVD.0<-function(OLS.list, Ob, r0=NULL){
  OLS.tsr=OLS.list$OLS.tsr
  p=OLS.tsr@modes[1]
  bc0=OLS.list$bc0
  summ.mat=OLS.list$summ.mat
  Theta=diag(0,p)
  tr.term<-diag(0,p)
  N<-0 #total sample size in the grps
  for(j in 1:nrow(Ob)){
    beta.j<-OLS.tsr@data[,Ob[j,1],Ob[j,2]]
    Theta=Theta+beta.j%*%t(beta.j)
    N=N+OLS.list$n.mat[Ob[j,1],Ob[j,2]]
  }
  Theta=Theta-bc0
  eigen.re<-eigen(Theta)
  if(is.null(r0)){
    thres<-sqrt(eigen.re$values[1]*max(p,log(N/nrow(Ob)))*nrow(Ob)/N)/2
    r0=sum(eigen.re$values>=thres)
  }
  eigen.re$vectors[,1:r0]
}

# ###DG-based transfer learning:aggregation with OLS
# TL.dg<-function(X.tar,y.tar, grp.ind, Beta.hat.tsr){
#   weight.hat<-NULL
#   for(ss in 1:30){
#     samp1=sample(1:length(y.tar),0.2*length(y.tar), replace=F)
#     ols.re=lm(y.tar[-samp1]~X.tar[-samp1,]-1)
#     ols.j=ols.re$coef
#     sig.hat=mean(ols.re$res^2)
#     dg.j=Beta.hat.tsr@data[,grp.ind[1],grp.ind[2]]
#     weight2=(sum((y.tar[samp1]-X.tar[samp1,]%*%dg.j)^2)-sum((y.tar[samp1]-X.tar[samp1,]%*%ols.j)^2))/2/sig.hat
#     weight.hat=c(weight.hat,1/(1+exp(weight2)))
#   }
#   weight=mean(weight.hat)
#   ols.j=lm(y.tar~X.tar-1)$coef
#   beta.tl=weight*dg.j+(1-weight)*ols.j
#   
#   beta.tl
# }
##DG-based transfer learning: Calibration with Trans-Lasso
TL.dg2<-function(X.tar,y.tar, grp.ind, Beta.hat.tsr){
  
  dg.j=Beta.hat.tsr@data[,grp.ind[1],grp.ind[2]]
  cv.re<-cv.glmnet(x=X.tar,y=y.tar-X.tar%*%dg.j)
  delta.hat<-coef(cv.re, s='lambda.min')[-1]
  beta.tl=dg.j+delta.hat
  
  beta.tl
}

TDplot<-function(B.hat,B.ora){
  p.vec<-B.hat@modes
  SSE.re<- NULL
  for(j in 1:p.vec[2]){
    for(k in 1:p.vec[3]){
      SSE.re<-rbind(SSE.re,c(j,k,
                             sum((B.hat@data[,j,k]-B.ora@data[,j,k])^2)))
    }
  }
  plot_ly(x=SSE.re[,1], y=SSE.re[,2], z=SSE.re[,3], type="scatter3d", mode="markers")
}




Maxmin<-function(beta.mat,Sigma){
  
  L=ncol(beta.mat)
  Gamma.positive <- diag(0,L)
  for(i in 1:L){
    for(j in 1:L){
      Gamma.positive[i,j]=t(beta.mat[,j])%*%Sigma%*%beta.mat[,i]
    }
  }
  eigen.val<-eigen(Gamma.positive)$values
  if(eigen.val[L]<0){
    Gamma.positive=Gamma.positive+diag(2*eigen.val[L],L)
  }
  v <- Variable(L)
  objective <- Minimize(quad_form(v, Gamma.positive))
  constraints <- list(v >= 0, sum(v) == 1)
  prob.weight <- Problem(objective, constraints)
  result <- solve(prob.weight)
  opt.weight= result$getValue(v)
  maxmin=beta.mat%*%opt.weight
  
  return(list(opt.weight=opt.weight, beta.mm=maxmin))
}
