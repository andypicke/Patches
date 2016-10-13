
##

# * to do: write a manipulate function to plto histogram of gamma with bars for adjusting thresholds for N2,dtdz,eps etc.

rm(list=ls())
the_file <- "/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/TestChiMethod/test.csv"
dat<-read.csv(the_file,header=FALSE,col.names = c("N2","dtdz","chi","eps"))

hist(log10(dat$N2))
hist(log10(dat$dtdz))
hist(log10(dat$chi))
hist(log10(dat$eps))

dat$gamma <- dat$N2 * dat$chi /2 / dat$eps / (dat$dtdz^2)
which(dat$gamma>100)
hist(log10(dat$gamma))

dat2<-subset(dat,log10(eps)>-8)
hist(log10(dat$gamma))
hist(log10(dat2$gamma))

library(ggplot2)
qplot(log10(dat$eps))

g<-ggplot(data=dat,aes(x=log10(dat$chi),y=log10(dat$eps)))
g+geom_density2d()
g+stat_bin2d()
g+stat_binhex()
g+stat_bin_hex()

smoothScatter(log10(dat$chi),log10(dat$eps),xlim=c(-13,-4))
