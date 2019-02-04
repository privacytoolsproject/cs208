##
##  hello.r
##
##  A simple demonstration of screen output, subsetting, and types in R
##
##  2019.2.4  jH
##


#### Printing ####

print("hello world")			# Print a line of text

							# End program at this line

h<-matrix(1,nrow=5,ncol=2)		# Create a matrix of numbers

print(h)						# Print matrix h

h<-"hello world"				# Overwrite h

print(h)						# Print h again





#### Creating Sequences ####

s<-1:15								 # Here is a simple way to construct a vector

print(s)

h<-matrix(s,nrow=5,ncol=3)      	 # We can reform the vector as a matrix

print(h)


t<-seq(from=1,by=3,length=5)    	 # Here is a more flexible function to construct a vector

flag<- round(t/2) == t/2			 # Flag is a logical (Boolean) vector with value TRUE if t is even

print(flag)

print(t[flag])                  	 # Very usefully you can subset using a boolean vector



print(h[2:4,3])						 # Matrices get subsetted first by row, then by column


#### Creating Loops ####

for(i in 1:5){                  	 # The indicator "i" iterates through the sequence 1,2,3,4,5
  print("The next element of t is")  # Everything inside the soft parentheses is repeated 
  print(t[i])                        #   every iteration with the current value of i
  thisExists <- TRUE				 # Note: Loops in R do not have their own scope
}

print(thisExists)




#### Simple Functions ####

something <- TRUE

myfunction<-function(x=1, y, ...){
	sum <- x+y
	something <- TRUE
	return(sum)
}


a<- 5
b<- 7
c<- 11

sum.1 <- myfunction(x=a, y=b)
print(sum.1)

sum.2 <- myfunction(y=7, z="notReallyAnArgument")
print(sum.2)




#### Dataframes ####

n<-20
a<- rep(c(TRUE,FALSE),n/2)
x<- rnorm(n=n, mean=0, sd=1)
y<- rbinom(n=n, size=1, prob=0.3)
z<- c("aa","bb","cccc","dddd")

mydata <- data.frame(a,x,y,z,b,c)

print(mydata)

stop()

# Methods of subsetting
# print(mydata$x)
# print(mydata[1:3,"x"])
# print(mydata[a,])
# print( names(mydata) )


#### Lists ####

mylist <- list( a, x, y)
catchphrase <- "differential privacy"
mylist[[5]] <- catchphrase


# print(list)
# print( list[[2]] )
# print( list$catchprase )
# print( names(mydata) )







