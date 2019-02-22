##
##  expEpsilon.r
##
##  demonstrate the approximation (1+\epsilon) \approx e^{\epsilon} for small \epsilon
##
##  jH 2019.2.22
##


epsilon<-seq(from=1, to=0.0001, length=20)
plot(1+epsilon, exp(epsilon), type="l", lwd=2, col="dodgerblue")
abline(a=0, b=1, lty=2, lwd=2)

dev.copy2pdf(file="./figs/expEpsilon.pdf")
