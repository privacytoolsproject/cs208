##
##  showDistributions.r
##
##  show probability distributions, cumulatives and their inverses, for Laplace distribution
##
##  jH 19.02.22
##

par(mfrow=c(3,2))
set.seed(123)

x<-seq(from=-3,to=3,length=100)

plot(x,dlap(x),type="l",main="density function")
plot(x,plap(x),type="l",main="cumulative density")

x<-seq(from=.001,to=.999,length=100)
plot(x,qlap(x),type="l",ylim=c(-3,3),main="inverse cumulative")

hist(rlap(size=1000),main="histogram of random draws", breaks=-9:9)
hist(rlap(size=1000),main="histogram of random Laplace draws", breaks=-9:9)
hist(qlap(runif(1000)),main="histogram of inv.cml.of random uniforms", breaks=-9:9)

dev.copy2pdf(file="./figs/laplaceDistributions.pdf")
