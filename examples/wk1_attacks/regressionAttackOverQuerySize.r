##
##  regressionAttackOverQuerySize.r
##
##  demonstrate reconstruction attack by regression on sums of subsets
##  cycle over different query sizes and show average correct reconstruction rate
##
##  jH 2019.2.5
##

#### Parameters ####
set.seed(123)
n <- 100        	# Dataset size
#k.trials <- 100  	# Number of queries
q.size <- 30    	# Subset size   
addPrior <- TRUE	# Whether to add priors
k.seq <- seq(from=40, to=150, by=1) # Sequence of query sizes
sims <- 50							# Number of iterations of each query size


#### Get Data ####
## If we wanted to simulate data, we might try:
#my.pi <- 0.2    
#sensitiveData <- rbinom(n, size=1, prob=my.pi)

## But here we read it from the repository:
pums <- read.csv(file="../../data/PUMS5extract10000.csv")
var <- "latino"    											# or also try "educ" for multivalued education scale
my.pi <- mean(pums[,var])
sampleIndex <- sample(x=1:nrow(pums), size=n, replace = FALSE )
sensitiveData <- pums[sampleIndex, var]


#### Here is our seemingly innocuous aggregated query ####
query <- function(n, data, sd=0){
	index <- sample(1:length(data), n)
	subset <- data[index]
	sum <- sum(subset) + rnorm(n=1, mean=0, sd=sd)     		# default is no noise, but can add
	return(list(sum=sum, index=index))
}


historyQuerySize <- matrix(NA, nrow=length(k.seq), ncol=3)

for(j in 1:length(k.seq)){       # Loop over number of queries
	if(round((j-1)/10) == ((j-1)/10)) cat(j, "of", length(k.seq), "\n")  # This is a longer simulation, so nice to have screen output showing progress
	k.trials <- k.seq[j]

	historyTemp <- matrix(NA, nrow=sims, ncol=2)

	for(k in 1:sims){

		#### Here we run our query repeatedly and record results ####
		history <- matrix(NA, nrow=k.trials, ncol=n+1)            	# a matrix to store simulation results in

		for(i in 1:k.trials){
			res <- query(n=q.size, data=sensitiveData, sd=0.25)
			indicator <- 1:n %in% res$index                         # convert indices into a series of boolean/dummy variables
			indicator <- as.numeric(indicator)
			history[i,] <- c(res$sum, indicator)                    # save into our results matrix
		}

		#### Add (data augmentation) prior ####

		if(addPrior){
			s <- max(100 - k.trials, 0)
			x <- matrix(rbinom(s*n, size=1, prob=0.5), nrow=s, ncol=n)
			y <- apply(x, 1, sum)*my.pi
			prior <- matrix(NA, nrow=s, ncol=n+1)
			prior[,1] <- y
			prior[,2:(n+1)] <- x
			history <- rbind(history, prior)
		}

		#### Convert matrix into data frame ####
		xnames <- paste("x", 1:n, sep="")
		varnames<- c("y", xnames)
		releaseData <- as.data.frame(history)                     # convert matrix into data frame
		names(releaseData) <- varnames                   


		#### Run a linear regression ####
		formula <- paste(xnames, collapse=" + ")                  # construct the formula, y ~ x1 ... xn -1
		formula <- paste("y ~ ", formula, "-1")
		formula <- as.formula(formula)

		output <- lm(formula, data=releaseData)                   # run the regression
		estimates <- output$coef                                  # save the estimates

		#### Determine which observations correctly reconstructed ####

		true.1 <- (estimates>0.5) & (sensitiveData==1)            # Correctly predicted values
		true.0 <- (estimates<0.5) & (sensitiveData==0)
		true1.frac <- round(sum(true.1, na.rm=TRUE)/sum(sensitiveData)*100)/100
		true0.frac <- round(sum(true.0, na.rm=TRUE)/sum(1-sensitiveData)*100)/100

		historyTemp[k,1] <- true1.frac   # Store the true fractions for this simulations
		historyTemp[k,2] <- true0.frac

	}

	historyQuerySize[j,1] <- k.seq[j]     # Store the average true fraction across simulations for this query size
	historyQuerySize[j,2] <- mean(historyTemp[,1])
	historyQuerySize[j,3] <- mean(historyTemp[,2])

}




#### Plot results ####
semi.blue <- rgb(0,90,239,200,maxColorValue=255)          # Slightly transparent colors
semi.red  <- rgb(239,90,0,200,maxColorValue=255)
col.values <- c(semi.red, semi.blue)

plot(x=historyQuerySize[,1], y=historyQuerySize[,2], ylim=c(0,1.2), type="l", lwd=1.5, col=semi.blue, main="Reconstruction Fraction Over Number of Queries", ylab="Correct Reconstruction Fraction", xlab="Number of Queries")
lines(x=historyQuerySize[,1], y=historyQuerySize[,3], type="l", lwd=1.5, col=semi.red)
abline(h=1, lty=2)

legx <- k.seq[round(length(k.seq)*.7)]

lines(x=c(legx, legx*.9), y=c(0.2,0.2), col=semi.red)
lines(x=c(legx, legx*.9), y=c(0.3,0.3), col=semi.blue)

text(x=legx, y=0.2, pos=4, labels="reconstructed 0's")
text(x=legx, y=0.3, pos=4, labels="reconstructed 1's")

#### Export graph to .pdf ####
dev.copy2pdf(file="./figs/regAttack_Over_Query_Size.pdf")





