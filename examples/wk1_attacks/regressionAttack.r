##
##  regressionAttack.r
##
##  demonstrate reconstruction attack by regression on sums of subsets
##
##  jH 2019.2.4
##

#### Parameters ####
set.seed(123)
n <- 100        # Dataset size
k.trials <- 80  # Number of queries
q.size <- 30    # subset size   
addPrior <- FALSE


#### Get Data ####
## If we wanted to simulate data, we might try:
#my.pi <- 0.2    
#sensitiveData <- rbinom(n, size=1, prob=my.pi)

## But here we read it from the repository:
pums <- read.csv(file="../../data/PUMS5extract10000.csv")
var <- "latino"
my.pi <- mean(pums[,var])
sampleIndex <- sample(x=1:nrow(pums), size=n, replace = FALSE )
sensitiveData <- pums[sampleIndex, var]


#### Here is our seemingly innocuous aggregated query ####
query <- function(n, data){
	index <- sample(1:length(data), n)
	subset <- data[index]
	sum <- sum(subset) #+ rnorm(n=1, mean=0, sd=0.1)
	return(list(sum=sum, index=index))
}


#### Here we run our query repeatedly and record results ####
history <- matrix(NA, nrow=k.trials, ncol=n+1)            # a matrix to store results in

for(i in 1:k.trials){
	res <- query(n=q.size, data=sensitiveData)
	indicator <- 1:n %in% res$index                         # convert indices into a series of boolean/dummy variables
	indicator <- as.numeric(indicator)
	history[i,] <- c(res$sum, indicator)                    # save into our results matrix
}

#### Add (data augmentation) prior
s <- max(100 - k.trials, 0)
x <- matrix(rbinom(s*n, size=1, prob=0.5), nrow=s, ncol=n)
y <- apply(x, 1, sum)*my.pi

if(addPrior){
	prior <- matrix(NA, nrow=s, ncol=n+1)
	prior[,1] <- y
	prior[,2:(n+1)] <- x
	history <- rbind(history, prior)
}

#### Convert matrix into data frame
xnames <- paste("x", 1:n, sep="")
varnames<- c("y", xnames)
releaseData <- as.data.frame(history)                     # convert matrix into data frame
names(releaseData) <- varnames                   


#### Run a linear regression ####
formula <- paste(xnames, collapse=" + ")                  # construct the formula, y ~ x1 ... xn -1
formula <- paste("y ~ ", formula, "-1")
formula <- as.formula(formula)
print(formula)

output <- lm(formula, data=releaseData)                   # run the regression
estimates <- output$coef                                  # save the estimates


#### Plot results ####
delta <- 0.05                                             # Slight disturbance to add
jitterx <- runif(n=n, min=-delta, max=delta)
jittery <- runif(n=n, min=-delta, max=delta)
semi.blue <- rgb(0,90,239,200,maxColorValue=255)          # Slightly transparent colors
semi.red  <- rgb(239,90,0,200,maxColorValue=255)
col.values <- c(semi.red, semi.blue)

true.1 <- (estimates>0.5) & (sensitiveData==1)            # Correctly predicted values
true.0 <- (estimates<0.5) & (sensitiveData==0)
true1.frac <- round(sum(true.1, na.rm=TRUE)/sum(sensitiveData)*100)/100
true0.frac <- round(sum(true.0, na.rm=TRUE)/sum(1-sensitiveData)*100)/100
truth.col<-1 + true.1 + true.0

plot(x=estimates + jitterx, y=sensitiveData + jittery, xlab="estimate", ylab="sensitive value", main="Reconstruction of Latino Variable", col=col.values[truth.col])    # Plot reconstruction against actual sensitive data
abline(v=0.5, lty=2)
text(x=0.5, y=0.8, labels=paste("fraction ones correct: ", true1.frac), pos=4)
text(x=0.5, y=0.2, labels=paste("fraction zeros correct: ", true0.frac), pos=2)

#### Export graph to .pdf
dev.copy2pdf(file="./figs/regAttack.pdf")





