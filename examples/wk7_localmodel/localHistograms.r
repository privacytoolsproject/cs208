##
##  localRelease.r
##
##  demonstrate local model of differential privacy
##
##  jH 2019.3.29
##

rm(list=ls())		# Remove any objects in memory
set.seed(123)

#### Release Functions ####

localRelease <- function(x, values=c(-1,1), epsilon){
	draw <- runif(n=1, min=0, max=1)
	cutoff <- 1/(1+exp(epsilon))
	if(draw<cutoff){
		return(values[!values%in%x])		
	}else{
		return(x)
	}
}


correction <- function(release, epsilon){
	inflation <- (exp(epsilon) + 1)/(exp(epsilon) - 1)
	expectation <- mean(release * inflation)
	return(expectation)
}


correction01 <- function(release, epsilon, sensitivity=1){
	inflation <- (exp(epsilon/sensitivity) + 1)/(exp(epsilon/sensitivity) - 1)
	release.trans <- (release-0.5)*2
	expectation <- release.trans * inflation
	expectation.trans <- expectation/2 + 0.5
	return(expectation.trans)
}



# Random draw from Laplace distribution
#
# mu numeric, center of the distribution
# b numeric, spread
# size integer, number of draws
# 
# return Random draws from Laplace distribution
# example:
# 
# rlap(size=1000)

rlap = function(mu=0, b=1, size=1) {
    p <- runif(size) - 0.5
    draws <- mu - b * sgn(p) * log(1 - 2 * abs(p))
    return(draws)
}

# Sign function
# 
# Function to determine what the sign of the passed values should be.
#
# x numeric, value or vector or values
# return The sign of passed values
# example:
#
# sgn(rnorm(10))

sgn <- function(x) {
    return(ifelse(x < 0, -1, 1))
}

## Bound/Censor/Clip a variable to a range
clip <- function(x, lower, upper){
	x.clipped <- x
	x.clipped[x.clipped<lower] <- lower
	x.clipped[x.clipped>upper] <- upper
	return(x.clipped)	
}

## Sample with replacement from a vector
bootstrap <- function(x, y=NULL, n){
	index <- sample(x=1:length(x), size=n, replace=TRUE) 

	if(is.null(y)){
		return(x[index])
	}else{
		return(list(x=x[index], y=y[index]))
	}
}



library("foreign")
PUMSdata <- read.csv(file="../../data/FultonPUMS5full.csv")   

data <- PUMSdata$educ    		# variable for means



## Differentially private histogram for integers
localHistogramRelease <- function(x, lower, upper, nbins=0, epsilon){
	n <- length(x)
	if(nbins==0){
		lower <- floor(lower)
		upper <- ceiling(upper)
		bins <- lower:upper   
        nbins <- length(bins)
    }

    x.clipped <- clip(x=x, lower=lower, upper=upper)

	sensitivity <- 2
	scale <- sensitivity / epsilon

	sensitiveValue <- DPrelease <- rep(NA,nbins)
	for(i in 1:length(bins)){
		sensitiveValue[i] <- sum(x.clipped==bins[i])
		DPrelease[i] <- localRelease(sensitiveValue[i], values=c(0,1), epsilon=epsilon/2)
	}


	return(list(release=DPrelease, true=sensitiveValue, codebook=bins))
}

nboot <- 20000
data1 <- bootstrap(data, n=nboot)






truefrac <- function(x, lower, upper){
	fractions <- hist(x, breaks=(lower:(upper+1)-0.5), plot=FALSE)$density
	return(fractions)
}

out1 <- matrix(NA, nrow=nboot, ncol=16)
for(i in 1:nboot){
	out1[i,] <- localHistogramRelease(x=data1[i], lower=1, upper=16, epsilon=0.5)$release
}


values <- apply(out1,2,mean)


par(mfcol=c(2,1))

barplot(values)
plot(values, type="h", lwd=3)





showHistLocal <- function(DPrelease, true, codebook, main="Histogram"){

	semi.blue <- rgb(0,90,239,150,maxColorValue=255)          # Slightly transparent colors
	semi.red  <- rgb(239,90,0,150,maxColorValue=255)

	allylim <- c(min(c(DPrelease,true), na.rm = TRUE), max(c(DPrelease, true), na.rm = TRUE))
	granularity <- (max(codebook) - min(codebook))/(length(codebook)-1)

	allxlim <- c(min(codebook) - 0.5*granularity, max(codebook + 0.5*granularity))

    # Build empty plot
	plot.new()
	plot.window( xlim=allxlim, ylim=allylim)
	title(main = main)
	axis( side=1 )
	axis( side=2 )

	tiny <- granularity*0.03 # slight spacing between bars
	overlap <- granularity*0.2 # some small overlap between sensitive and DP values

	for(i in 1:length(codebook)){
		rect(xleft=codebook[i]-overlap, ybottom=0, xright=codebook[i]+0.5*granularity-tiny, ytop=true[i], col=semi.red)
		rect(xleft=codebook[i]-0.5*granularity+tiny, ybottom=0, xright=codebook[i]+overlap, ytop=DPrelease[i], col=semi.blue)
	}
}


values <- apply(out1,2,mean)
DPmeans <- correction01(values, epsilon=0.5)
true <- truefrac(data1, lower=1, upper=16)

par(mfcol=c(1,1))
showHistLocal(DPrelease=DPmeans, true=true, codebook=1:16, main="Histogram of Local Model Release of Education")

dev.copy2pdf(file="./figs/localEducHist.pdf")

#### Actual libraries for hash functions
#install.packages("openssl")
library("openssl")
sha256(c("james","salil","james"), key="my_secret")


#### Terrible Hash Functions
# use first letter of string
thash <- function(x){
	x <- tolower(x)
	hash<-NULL
	for(i in 1:length(x)){
		first.letter <- substr(x[i], start=1, stop=1)
		temp <- which(first.letter==letters)      
		hash <- c(hash, max(temp,0))         # max helps map nonletters to 0
	}
	return(hash)
}

# use last letter of string
thash2 <- function(x){
	x <- tolower(x)
	hash<-NULL
	for(i in 1:length(x)){
		last.letter <- substr(x[i], start=nchar(x[i]), stop=nchar(x[i]))
		temp <- which(last.letter==letters)      
		hash <- c(hash, max(temp,0))         # max helps map nonletters to 0
	}
	return(hash)
}


#### Discover Species Names

data("iris")

names <- bootstrap(iris$Species, n=5000) 
names.hash <- thash(names)

out2 <- out3 <- matrix(NA, nrow=length(names), ncol=27)

for(i in 1:length(names.hash)){
	out2[i,] <- localHistogramRelease(x=names.hash[i], lower=0, upper=26, epsilon=0.5)$release
}


names.hash2 <- thash2(names)
for(i in 1:length(names.hash)){
	out3[i,] <- localHistogramRelease(x=names.hash2[i], lower=0, upper=26, epsilon=0.5)$release
}


values <- apply(out2,2,mean)
DPmeans <- correction01(values, epsilon=0.5)
true <- truefrac(names.hash, lower=0, upper=26)


showHistLocal(DPrelease=DPmeans, true=true, codebook=0:27, main="Histogram of Local Model Release of Name Hash")



#### Show Client-SFP for string discovery ####

clientSFP <- function(x, epsilon, myhash){
	a <- 1 # Correct this
	l <- 1 # Correct this
	b <- substr(x,start=l,stop=l)
	return(list(a=a, b=b, l=l))
}

x<- bootstrap(iris$Species, n=10000)
l <- rep(1,length(x))
b <- rep("a",length(x))
out4 <- matrix(NA, nrow=length(x), ncol=27)
myepsilon <- 2


for(i in 1:length(x)){
	release <- clientSFP(x[i], epsilon=myepsilon, myhash=thash)
	out4[i,] <- release$a
	b[i] <- release$b
	l[i] <- release$l
}


# Show identified hash

codebook<- 0:26
values <- apply(out4,2,mean)
DPmeans <- correction01(values, epsilon=myepsilon)
true <- truefrac(thash(x), lower=0, upper=26)

showHistLocal(DPrelease=DPmeans, true=true, codebook=0:27, main="Histogram of Local Model Release of Name Hash")
dev.copy2pdf(file="./figs/localHashHist.pdf")



# Piece together puzzle

Threshold <- 0.1
discovered <- which(DPmeans>Threshold) 

cat("Actual Names: \n")
cat(paste(sort(unique(iris$Species)), "\n"))

for(j in 1:length(discovered)){
	flag <- out4[,discovered[j]] == 1
	temp.b <- b[flag]
	temp.l <- l[flag]
	t <- table(temp.b,temp.l)
	print(t)
	size <- ncol(t)
	word <- rep("",size)
	for(k in 1:ncol(t)){
		word[k]<-row.names(t)[which(t[,k]==max(t[,k]))][1]
	}
	print(word)
}





