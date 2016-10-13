
# Using manipulate package to investigate effect of eps,n2,dtdz thresholds on distribution of gamma in eq14 data

rm(list=ls())
library(manipulate)

the_file <- "/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/TestChiMethod/test.csv"
dat<-read.csv(the_file,header=FALSE,col.names = c("N2","dtdz","chi","eps"))

hist(log10(dat$N2))
hist(log10(dat$dtdz))
hist(log10(dat$chi))
hist(log10(dat$eps))

dat$gamma <- dat$N2 * dat$chi /2 / dat$eps / (dat$dtdz^2)


# Histogram of gamma, sliders for eps min
g<-function(epsmin){
        dat2<-subset(dat,log10(eps)>epsmin)
        dat2<-subset(dat2,gamma>0 & gamma<2)
        hist(dat2$gamma,breaks=seq(0,2,0.025),xlim = c(0,0.4))
        abline(v=mean(dat2$gamma),col="red")
        abline(v=median(dat2$gamma),col="blue")
}

manipulate(g(epsmin),epsmin=slider(-10,-6,step=0.5))

# Histogram of gamma, sliders for eps,n2,min
g2<-function(epsmin,n2min){
        dat2<-subset(dat,log10(eps)>epsmin & log10(N2)>n2min)
        dat2<-subset(dat2,gamma>0 & gamma<2)
        hist(dat2$gamma,breaks=seq(0,2,0.025),xlim = c(0,0.4))
        abline(v=mean(dat2$gamma),col="red")
        abline(v=median(dat2$gamma),col="blue")
}

manipulate(g2(epsmin,n2min),epsmin=slider(-10,-6,step=0.5),n2min=slider(-6,-1))

# Histogram of gamma, sliders for eps,n2,dtdz min
g3<-function(epsmin,n2min,dtdzmin){
        dat2<-subset(dat,log10(eps)>epsmin & log10(N2)>n2min & log10(dtdz)>dtdzmin)
        dat2<-subset(dat2,gamma>0 & gamma<2)
        hist(dat2$gamma,breaks=seq(0,2,0.025),xlim = c(0,0.3))
        abline(v=mean(dat2$gamma),col="red")
        abline(v=median(dat2$gamma),col="blue")
}

manipulate(g3(epsmin,n2min,dtdzmin),epsmin=slider(-10,-6,step=0.5),n2min=slider(-6,-1,step=0.5),dtdzmin=slider(-6,-0,step=0.5))


# Make boxplot of gamma
g4<-function(epsmin,n2min,dtdzmin){
        dat2<-subset(dat,log10(eps)>epsmin & log10(N2)>n2min & log10(dtdz)>dtdzmin)
        dat2<-subset(dat2,gamma>0 & gamma<2)
boxplot(dat2$gamma)
}

manipulate(g4(epsmin,n2min,dtdzmin),epsmin=slider(-10,-6,step=0.5),n2min=slider(-6,-1,step=0.5),dtdzmin=slider(-6,-0,step=0.5))
