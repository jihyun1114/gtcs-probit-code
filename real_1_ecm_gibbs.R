# =============================================================================
# real_1_ecm_gibbs.R  — FINAL real-data analysis for Table 5 + Figure 3.
# CONFIRMED configuration:
#   * 2-covariate model: Gender, C-Amer   (NO interaction)   [advisor sign-off pending]
#   * probit, alpha=0.947, beta=0.989, random pool-of-5 (seed 2025)
#   * baseline: I-splines, m=10 quantile knots (ties nudged), grid age 0..100
#   * Gibbs baseline-intercept prior: gam0 ~ N(m0=-3, 1/v0=1/10)   [gives S(0)~=1]
#   * ECM re-run on the SAME pooled data for a consistent Table 5
# Outputs (realdata_results/): table5_ecm.csv, table5_gibbs.csv,
#   fig3_ecm.csv, fig3_gibbs.csv, realplot_final.png
# =============================================================================
source("ECM.R"); source("Gibbs.R")

ALPHA<-0.947; BETA<-0.989; POOL<-5; SEED<-2025; GRID<-seq(0,100,length.out=200); M<-10
MAXIT<-5000; BURN<-1000; M0<--3; V0<-10
OUT<-"realdata_results"; dir.create(OUT,showWarnings=FALSE)

raw<-read.table("g.2008.swab.xls",header=TRUE); keep<-raw$age>0
age<-raw$age[keep]; gender<-as.integer(raw$gender[keep]); camer<-as.integer(raw$race[keep]==1)
X<-cbind(Gender=gender, CAmer=camer)                         # 2 covariates, no interaction
C<-age; ind<-as.integer(raw$test_result[keep])
n<-nrow(X); set.seed(SEED); ord<-sample.int(n); X<-X[ord,];C<-C[ord];ind<-ind[ord]
ng<-floor(n/POOL); rem<-n-POOL*ng
groupID<-c(rep(seq_len(ng-1),each=POOL),rep(ng,POOL+rem)); groupsize<-c(rep(POOL,ng-1),POOL+rem)
Dtrue<-as.integer(tapply(ind,groupID,function(z)any(z==1)))
Ypool<-Dtrue; Ypool[Dtrue==1]<-rbinom(sum(Dtrue==1),1,ALPHA); Ypool[Dtrue==0]<-rbinom(sum(Dtrue==0),1,1-BETA)
Y<-rep(Ypool,times=groupsize)
KN<-as.numeric(c(0,quantile(C,seq(0.1,0.9,length.out=M)),max(C)+0.01))
for(i in 2:length(KN)) if(KN[i]<=KN[i-1]) KN[i]<-KN[i-1]+1e-2

cat("Running ECM ...\n")
fe<-ECM(Y=Y,X=X,c=C,groupID=groupID,alpha=ALPHA,beta=BETA,gam0=-3,
        gam=c(1,.7,.2,.1,.1,.1,rep(0,7)),grid=GRID,knots=KN)
cat("Running Gibbs (gam0 prior m0=-3, v0=10) ...\n")
set.seed(1)
fg<-Gibbs(Y=Y,X=X,c=C,groupID=groupID,alpha=ALPHA,beta=BETA,grid=GRID,knots=KN,
          gam0=-3,m0=M0,v0=V0,maxiter=MAXIT,burn.in=BURN)

nm<-c("Gender","C-Amer")
ecm.tab<-data.frame(Covariate=nm, Est=round(fe$summary.theta$EST,3),
  L95=round(fe$summary.theta$lb95,3), U95=round(fe$summary.theta$ub95,3))
th<-fg$theta.mat[-(1:BURN),]
gib.tab<-data.frame(Covariate=nm, Est=round(colMeans(th[,1:2]),3),
  L95=round(apply(th[,1:2],2,quantile,.025),3), U95=round(apply(th[,1:2],2,quantile,.975),3))
write.csv(ecm.tab,file.path(OUT,"table5_ecm.csv"),row.names=FALSE)
write.csv(gib.tab,file.path(OUT,"table5_gibbs.csv"),row.names=FALSE)
cat("\n===== TABLE 5 (FINAL: 2-cov, Gibbs gam0~N(-3,1/10)) =====\n--- ECM ---\n"); print(ecm.tab,row.names=FALSE)
cat("--- Gibbs ---\n"); print(gib.tab,row.names=FALSE)
cat(sprintf("Gibbs gam0=%.3f -> S(0)=%.3f    ECM gam0=%.3f -> S(0)=%.3f\n",
    mean(th[,3]),1-pnorm(mean(th[,3])), fe$gam0, 1-pnorm(fe$gam0)))

## Figure 3
grp<-list("Female, Other"=c(0,0),"Male, Other"=c(1,0),"Female, C-Amer"=c(0,1),"Male, C-Amer"=c(1,1))
g.ecm<-qnorm(fe$St,lower.tail=FALSE); S.ecm<-sapply(grp,function(x)1-pnorm(g.ecm+drop(x%*%fe$theta.est)))
g.gib<-colMeans(fg$g.mat[-(1:BURN),]); tmean<-colMeans(th[,1:2]); S.gib<-sapply(grp,function(x)1-pnorm(g.gib+drop(x%*%tmean)))
write.csv(data.frame(t=GRID,S.ecm,check.names=FALSE),file.path(OUT,"fig3_ecm.csv"),row.names=FALSE)
write.csv(data.frame(t=GRID,S.gib,check.names=FALSE),file.path(OUT,"fig3_gibbs.csv"),row.names=FALSE)
png(file.path(OUT,"realplot_final.png"),width=1000,height=500)
par(mfrow=c(1,2),mar=c(4,4,2,0.5)); sel<-GRID<80
matplot(GRID[sel],S.ecm[sel,],type="l",lwd=2,lty=1:4,col=1:4,ylim=c(0.8,1),xlab="Age",ylab=expression(S(t)),main="ECM")
legend("bottomleft",names(grp),lwd=2,lty=1:4,col=1:4,bty="n")
matplot(GRID[sel],S.gib[sel,],type="l",lwd=2,lty=1:4,col=1:4,ylim=c(0.8,1),xlab="Age",ylab=expression(S(t)),main="Gibbs")
legend("bottomleft",names(grp),lwd=2,lty=1:4,col=1:4,bty="n")
dev.off(); cat("saved table5_*.csv, fig3_*.csv, realplot_final.png\n")
