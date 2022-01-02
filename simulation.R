set.seed(0)

p<-0.5
X1<-matrix(rbinom(2000,2,p),nrow=1000)
X2<-matrix(rbinom(2000,10,p),nrow=1000)
X3<-matrix(rnorm(2000,0,1),nrow=1000)

X<-rbind(kronecker(rep(1,4),cbind(X1,X2,X3)),kronecker(rep(1,4),cbind(X1,X2,X3)))

Z<-kronecker(diag(2),kronecker(rep(1,4),diag(1000)))

betaTG<-rnorm(6,0,1)
betaHDL<-rnorm(6,0,1)

uTG<-rnorm(1000,0,1)
uHDL<-rnorm(1000,0,1)

rho1TG<-0.9
rho2TG<-0.9
rho3TG<-0.9

rho1HDL<-0.9
rho2HDL<-0.9
rho3HDL<-0.9

K<-matrix(1,dim<-c(1000,1000))
diag(K)<-1

sigmaTG<-1
sigmaHDL<-1

sigmauTG<-1
sigmauHDL<-1

rhouTGHDL<-0.5

rhoeTGHDL<--0.5

SigmaTG<-sigmaTG*matrix(c(1,rho1TG,rho1TG*rho2TG,rho1TG*rho2TG*rho3TG,rho1TG,1,rho2TG,rho2TG*rho3TG,rho1TG*rho2TG,rho2TG,1,rho3TG,rho1TG*rho2TG*rho3TG,rho2TG*rho3TG,rho3TG,1),nrow=4)

SigmaHDL<-sigmaHDL*matrix(c(1,rho1HDL,rho1HDL*rho2HDL,rho1HDL*rho2HDL*rho3HDL,rho1HDL,1,rho2HDL,rho2HDL*rho3HDL,rho1HDL*rho2HDL,rho2HDL,1,rho3HDL,rho1HDL*rho2HDL*rho3HDL,rho2HDL*rho3HDL,rho3HDL,1),nrow=4)

eiTG<-










