

## Here is the likelihood function for a Logit
calcllik<-function(b,data){           
  y<-data[,1]
  x<-data[,2]

  pi<- 1/(1+exp(-b[1] - b[2]*x))        # Here is the systematic component

  llik<-y * log(pi) + (1-y) * log(1-pi) # Here is the stocastic component
  return(-llik)
}

## Differentially private mean release
gaussianReleaseNoise <- function(size=1, sensitivity, epsilon, delta){
	scale <- sensitivity *log(1.25/delta)/ epsilon
	noise <- rnorm(n=size, mean=0, sd=scale)
	return(noise)
}

## Bound/Censor/Clip a variable to a range
clip <- function(x, lower, upper){
	x.clipped <- x
	x.clipped[x.clipped<lower] <- lower
	x.clipped[x.clipped>upper] <- upper
	return(x.clipped)	
}




#### Run with actual data

library("foreign")
PUMSdata <- read.csv(file="../../data/MaPUMS5full.csv")

mydata<-PUMSdata[c("married","educ")]

output <- glm(married ~ educ, family="binomial", data=mydata)

print(summary(output))


#### Show the estimated model

xseq <- seq(from=-40, to=60, length=100)
f <- 1/(1 + exp(-output$coef[1] -output$coef[2]*xseq))

par(mfcol=c(2,1))

plot(xseq, f, type="l", lwd=1.5, col="red", ylim=c(0,1), ylab="E(y|x,theta)", xlab="education", main="Probability Married by Education")
abline(v=1, col="blue", lty=3)
abline(v=16, col="blue", lty=3)

plot(xseq, f, type="l", lwd=1.5, col="red", ylim=c(0,1), ylab="E(y|x,theta)", xlab="education", xlim=c(-5,20))
for(i in 1:16){
	flag<-mydata$educ==i
	points(x=i, y=mean(mydata[flag,"married"]))
}
dev.copy2pdf(file="./figs/married_educ.pdf")


#### Show the LogLikelihood surface

sample.data <- mydata[sample(1:nrow(mydata),10000), ]
b1.seq <- seq(from=-3, to=2, length=25)
b2.seq <- seq(from=-.2, to=.3, length=25)
llsurface <- matrix(NA, nrow=length(b1.seq), ncol=length(b2.seq))

for(i in 1:length(b1.seq)){
	for(j in 1:length(b2.seq)){
		llsurface[i,j] <- sum(-1* calcllik(b=c(b1.seq[i], b2.seq[j]), data=sample.data) )
	}
}

filled.contour(x=b1.seq, y=b2.seq, z=llsurface, color = terrain.colors,  xlab="constant parameter", ylab="education parameter")

dev.copy2pdf(file="./figs/logitLLike.pdf")




calcgradient <- function(B, C, theta, fun){
	dx <- 0.0001

	out1 <-	eval(fun(b=theta, data=B))
	out2 <- eval(fun(b=theta + c(0,dx), data=B))
	out3 <- eval(fun(b=theta + c(dx,0), data=B))

	Del.1 <- 1
	# Del.1 <- clip(Del.1, lower=0, upper=1)  # Fix this
	mean.Del.1 <- mean(Del.1)

	Del.2 <- 1
	# Del.2 <- clip(Del.2, lower=0, upper=1)  # Fix this
	mean.Del.2 <- mean(Del.2)

	return(c(mean.Del.1,mean.Del.2))
}



N <- nrow(mydata)
L <- round(sqrt(nrow(mydata)))

steps <- 10   	  # Fix this

## Shuffle the data
index <- sample(1:nrow(mydata))
mydata <- mydata[index,]
epsilon <-1

theta <- c(0,0)   # Starting parameters
C <- 10			  # Interval to clip over
nu <- c(1,0.01)   # Learning speeds


history <- matrix(NA, nrow=steps+1, ncol=2)
history[1,] <- theta

for(i in 1:steps){
	startB <- ((i-1)*L+1)
	if(i<L){
		stopB <- i*L
	}else{
		stopB <- nrow(mydata)
	}

	index<-sample(1:nrow(mydata),L)
	B <- mydata[startB:stopB, ]
	Del <- calcgradient(B, C, theta, fun=calcllik)
	cat("Del:  ",Del,"\n")
	theta <- theta   				# Fix this
	cat("Theta:",theta, "\n")

	history[i+1,] <- theta

}

par(mfcol=c(2,1))

all.ylim<-c( min(c(history[,1],output$coef[1] )), max(c(history[,1],output$coef[1] )))
plot(history[,1], type="l", ylim=all.ylim, ylab="beta 0", xlab="step", lwd=1.5)
abline(h=output$coef[1], lty=2, col="blue", lwd=1.5)


all.ylim<-c( min(c(history[,2],output$coef[2] )), max(c(history[,2],output$coef[2] )))
plot(history[,2], type="l", ylim=all.ylim, ylab="beta 1", xlab="step", lwd=1.5)
abline(h=output$coef[2], lty=2, col="blue", lwd=1.5)

dev.copy2pdf(file="./figs/dpSGD.pdf")



