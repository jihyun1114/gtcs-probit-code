# Fang final: fit + profile-likelihood SE with the I-spline matched to ECM/Gibbs (order=3, n.int=10).
source("fang_code/Support_functions_fixed.R")
suppressMessages({library(doParallel); library(foreach)})

ALPHA<-0.947; BETA<-0.989; POOL<-5; SEED<-2025; ORD<-3; NINT<-10
raw<-read.table("g.2008.swab.xls",header=TRUE); keep<-raw$age>0
age<-raw$age[keep]; gender<-as.integer(raw$gender[keep]); camer<-as.integer(raw$race[keep]==1)
z<-cbind(gender,camer); xt<-age; ind<-as.integer(raw$test_result[keep])
n<-nrow(z); set.seed(SEED); ord<-sample.int(n); z<-z[ord,]; xt<-xt[ord]; ind<-ind[ord]
ng<-floor(n/POOL); rem<-n-POOL*ng
groups<-c(rep(seq_len(ng-1),each=POOL),rep(ng,POOL+rem)); g.n<-c(rep(POOL,ng-1),POOL+rem)
upool.id<-unique(groups)
Dtrue<-as.integer(tapply(ind,groups,function(zz)any(zz==1)))
Ypool<-Dtrue; Ypool[Dtrue==1]<-rbinom(sum(Dtrue==1),1,ALPHA); Ypool[Dtrue==0]<-rbinom(sum(Dtrue==0),1,1-BETA)
Y<-Ypool; nu<-rep(ALPHA,length(Y)); omega<-rep(BETA,length(Y)); p<-ncol(z)
registerDoParallel(cores=p)

t0<-proc.time()
fit<-probit_GT(xt=xt,groups=groups,g.n=g.n,n.int=NINT,order=ORD,z=z,upool.id=upool.id,Ys=Y,nu=nu,omega=omega)
h<-5*length(upool.id)^(-0.5); E<-h*diag(p)
Me<-foreach(i=1:p,.combine=rbind) %dopar% {
  pl1<-prolik(xt=xt,groups=groups,g.n=g.n,n.int=NINT,order=ORD,z=z,upool.id=upool.id,Ys=Y,be=fit$beta+E[,i],g.est=fit$g,nu=nu,omega=omega)
  pl2<-prolik(xt=xt,groups=groups,g.n=g.n,n.int=NINT,order=ORD,z=z,upool.id=upool.id,Ys=Y,be=fit$beta-E[,i],g.est=fit$g,nu=nu,omega=omega)
  (pl1-pl2)/(2*h)
}
I<-Me%*%t(Me); se<-sqrt(diag(solve(I))); beta<-as.vector(fit$beta)
tab<-data.frame(Covariate=c("Gender","C-Amer"), Est=round(beta,3), SE=round(se,3),
                L95=round(beta-1.96*se,3), U95=round(beta+1.96*se,3),
                p=round(2*(1-pnorm(abs(beta)/se)),3))
cat(sprintf("Fang (order=%d, n.int=%d), %.0fs\n", ORD, NINT, (proc.time()-t0)["elapsed"]))
print(tab, row.names=FALSE)
write.csv(tab, "realdata_results/table5_fang.csv", row.names=FALSE)
