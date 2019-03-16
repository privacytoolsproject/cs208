# Install dependency packages if not current

james_update_packages <- function(packageList){
	availableRepos <- getCRANmirrors()
	flag <- availableRepos$Country=="USA" & grepl("https",availableRepos$URL,)
	useRepos <- sample(availableRepos$URL[flag],1)

	## install missing packages, and update if newer version available
	for(i in 1:length(packageList)){
		if (!require(packageList[i],character.only = TRUE)){
    		install.packages(packageList[i], repos=useRepos)
   		}
	}

	update.packages(ask = FALSE, dependencies = c('Suggests'), oldPkgs=packageList, repos=useRepos)
}

packagelist <- c("devtools", "jsonlite", "openssl")
james_update_packages(packagelist)


# Install PSIlence from GitHub
devtools::install_github("privacytoolsproject/PSI-Library", ref="develop") 
library("PSIlence")


# Loading test data
data(PUMS5extract10000) # Load test dataset of PUMS responses, 10000 samples from California

my_data <- PUMS5extract10000 # This is included in the PSI library

my_data$race <- "white"   # Add categorical variable constructed from binary indicators
my_data$race[my_data$black==1] <- "black"
my_data$race[my_data$asian==1] <- "asian"
my_data$race[my_data$latino==1] <- "latino"

# Some necessary metadata
racevalues <- c("white", "black", "asian", "latino")
my_n <- nrow(my_data)  # Dataset size
my_age_bounds <- c(0,110)
my_income_bounds <- c(0,100000)

## Generate and print dp release of counts of each type
dp.release1 <- dpHistogram$new(mechanism="mechanismLaplace", var.type="character", variable="race", n=my_n, epsilon=0.5, bins=racevalues, n.bins=length(racevalues))
dp.release1$release(my_data)
print(dp.release1$result$release)


## Generate and print dp release of mean of age
dp.release2 <- dpMean$new(mechanism='mechanismLaplace', var.type='numeric', variable="age", n=my_n, epsilon=0.1,  rng=my_age_bounds)
dp.release2$release(my_data)

mean.CI <- mean.getCI(release=dp.release2$result$release, epsilon=dp.release2$epsilon, sensitivity=diff(dp.release2$rng)/dp.release2$n, alpha=0.05)
mean.epsilon <- mean.getParameters(accuracy=0.1, n=my_n, alpha=0.05, rng=my_age_bounds)
mean.errorPromise <- mean.getAccuracy(epsilon=0.1, n=my_n, alpha=0.05, rng=my_age_bounds)
mean.JSON <- PSIlence:::mean.getJSON()

cat("release:\n", dp.release2$result$release, "\n")
cat("CI:\n", mean.CI, "\n")
cat("error promise:\n", mean.epsilon, "\n")
cat("epsilon:\n", mean.errorPromise, "\n")
cat("empty JSON:", mean.JSON, "\n")


## Generate and print binary tree of income
dp.release3 <- dpTree$new(mechanism="mechanismLaplace", var.type="numeric", variable="income", n=my_n, rng=my_income_bounds, gran=100000/32, epsilon=0.3, alpha=0.05)
dp.release3$release(my_data)
print(dp.release3)


## Show use of optimal composition theorem from above releases

globalDelta <- 10^-6
eps.1 <- 0.5
eps.2 <- 0.1
eps.3 <- 0.3
delta.1 <- delta.2 <- delta.3 <- globalDelta/6   # Note we assigning are using delta_i = globalDelta/(2*k) for each of our k queries, 
												  # and leaving the rest of the global delta available to the optimal composition function

params <- matrix(c(eps.1, delta.1,
	               eps.2, delta.2,
	               eps.3, delta.3), nrow=3, ncol=2, byrow=TRUE)

# The following function exists in the package to give composition by optimal composition theorem
# ':::'' is used for a function that is not exported from the package
# For implementation see: https://github.com/privacytoolsproject/PSI-Library/blob/develop/R/CompositionTheorems.R
# Or from R:> print(PSIlence:::KOVhet)
# 
# Args:
#	params: a kx2 matrix where the first column corresponds to epsilon_i values and the second 
# 			corresponds to delta_i values. 
#	d_g: global delta value
#   print: Boolean that if TRUE will print each individual term of the theorem rather than just
#          the minimimum.
#   
# Returns:
#	global epsilon value guaranteed from the composition

out <- PSIlence:::KOVhet(params=params, d_g=globalDelta, print=TRUE)  


# And the following function gives the inverse, that is, given global epsilon and delta, finds a best (here equal) division
# among k releases that under optimal composition satisfies the global totals
# ':::'' is used for a function that is not exported from the package
# For implementation see: https://github.com/privacytoolsproject/PSI-Library/blob/develop/R/update_parameters.R
# Or from R:> print(PSIlence:::update_parameters)
#
# Args:
#	params: kx2 matrix of privacy parameters where column one corresponds
#			to epsilons and column two is deltas.
#	hold: vector of indices corresponding to rows of params that will not 
#		   be updated, either because they were just added or because the 
#		   user has requested that these values stay fixed (Hold feature). 
#	       If we are to update every parameter, set hold to 0. 
#	eps: global epsilon
#	del: global delta
#
# Returns:
#	kx2 matrix of updated parameters

k <- 3
epsilonGlobal <- 1
deltaGlobal <- 1e-9
init <- rep(c(1/k, 0), k )

params <- matrix(init, nrow=k, ncol=2, byrow=TRUE)

inverse <- PSIlence:::update_parameters(params=params, hold=0, eps=epsilonGlobal, del=deltaGlobal)

print(inverse)


## Demonstrate Error Promises and Epsilon Calculations

# Metadata values:
my.seq <- seq(from=log10(200), to=log10(1500), length=20)  	# make evenly spaced in logarithmic space
n.seq  <- round(10^my.seq)                                 	# round to integers

my.seq <- seq(from=log10(1), to=log10(0.01), length=5)     	# make evenly spaced in logarithmic space
ep.seq <- round(10^my.seq * 100) /100						# round to two decimal places

my.seq <- seq(from=log10(0.1), to=log10(0.01), length=5)
acc.seq <- round(10^my.seq * 100) /100

myrng <- c(0,1)

# Storage matrix:
agghistory <- matrix(NA, nrow=length(n.seq)*length(ep.seq), ncol=5)         # matrix to store results
aggcount <- 0                                               # counter

# Simulation:
for(i in 1:length(n.seq)){
	for(j in 1:length(ep.seq)){
		aggcount <- aggcount + 1
		agghistory[aggcount,1] <- n.seq[i]
		agghistory[aggcount,2] <- ep.seq[j]
		agghistory[aggcount,3] <- mean.getAccuracy(epsilon=ep.seq[j], n=n.seq[i], alpha=0.05, rng=myrng)
		agghistory[aggcount,4] <- mean.getParameters(accuracy=acc.seq[j], n=n.seq[i], alpha=0.05, rng=myrng)

	}
}


## Graphs for plotting results

par(mfrow=c(2,2))
color.palette<-rainbow(length(ep.seq), start=.7, end=.1)   # This creates a sequence of colors to use in subsequent plots, as in showchisq.r


for(j in 1:length(ep.seq)){
	flag <- agghistory[,2] == ep.seq[j]
	subhistory <- agghistory[flag,]

	allylim <- c(0, max(agghistory[,3]))

	if(j==1){
		plot(subhistory[,1],subhistory[,3], ylim=allylim, type="l", col=color.palette[j], xlab="N", ylab="Error Promise")
	}else{
		lines(subhistory[,1],subhistory[,3], col=color.palette[j])
	}
}

for(j in 1:length(ep.seq)){
	flag <- agghistory[,2] == ep.seq[j]
	subhistory <- agghistory[flag,]

	allylim <- c(min(agghistory[,3]), max(agghistory[,3]))

	xloc <- round(length(n.seq)*0.3)

	if(j==1){
		plot(subhistory[,1],subhistory[,3], ylim=allylim, type="l", log = "y", col=color.palette[j], xlab="N", ylab="Error Promise")
		text(x=subhistory[xloc,1], y=subhistory[xloc,3], label=  bquote(paste(epsilon == .(ep.seq[j]))), col=color.palette[j], pos=4)
	}else{
		lines(subhistory[,1],subhistory[,3], col=color.palette[j])
		text(x=subhistory[xloc,1], y=subhistory[xloc,3], label=  bquote(paste(epsilon == .(ep.seq[j]))), col=color.palette[j], pos=4) 
	}

}

for(j in 1:length(ep.seq)){
	flag <- agghistory[,2] == ep.seq[j]
	subhistory <- agghistory[flag,]

	allylim <- c(0, max(agghistory[,4]))

	if(j==1){
		plot(subhistory[,1],subhistory[,4], ylim=allylim, type="l", col=color.palette[j], xlab="N", ylab="Epsilon")
	}else{
		lines(subhistory[,1],subhistory[,4], col=color.palette[j])
	}
}

for(j in 1:length(ep.seq)){
	flag <- agghistory[,2] == ep.seq[j]
	subhistory <- agghistory[flag,]

	allylim <- c(min(agghistory[,4]), max(agghistory[,4]))

	xloc <- round(length(n.seq)*0.3)

	if(j==1){
		plot(subhistory[,1],subhistory[,4], ylim=allylim, type="l", log = "y", col=color.palette[j], xlab="N", ylab="Epsilon")
		text(x=subhistory[xloc,1], y=subhistory[xloc,4], label=  bquote(paste("Error" == .(acc.seq[j]))), col=color.palette[j], pos=4)
	}else{
		lines(subhistory[,1],subhistory[,4], col=color.palette[j])
		text(x=subhistory[xloc,1], y=subhistory[xloc,4], label=  bquote(paste("Error" == .(acc.seq[j]))), col=color.palette[j], pos=4) 
	}

}


dev.copy2pdf(file="./figs/psiPromises.pdf")


