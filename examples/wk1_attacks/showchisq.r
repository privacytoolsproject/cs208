##
##  showchisq.r
##
##  show the Chi^2 distribution function, for differing degrees of freedom
##  jH, 1/18/08
##


x<-1:1000/100                                     # This is going to make a finely graded vector to calculate through

maxdf<-5                                          # The largest value of the df parameter we want to graph

color.palette<-rainbow(maxdf, start=.7, end=.1)   # This creates a sequence of colors to use in subsequent plots

for(i in 1:maxdf){                                                             # Cycle through df values
  plot(x,dchisq(x,df=i),ylim=(c(0,1)),col=color.palette[i],type="l",lwd=2 )    # Plot the density (d-)
  abline(v=qchisq(0.9,df=i),col=color.palette[i],lty=2)                        # Overlay the inv. cumulative (q-)
  par(new=TRUE)                                                                # Keep the graph open for further plots
}
par(new=FALSE)								       # Close the graph
