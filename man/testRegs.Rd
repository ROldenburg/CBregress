\name{testRegs}
\title{Test regression methods on sampled data}
\description{
This function samples data and applies standard, orthogonal and case based linear regression. 
Standard deviation of the error between estimated and true values for the slopes are given.}
\usage{
testRegs(m,n,A,sig1,sig2,generator=runif)
}
\arguments{
\item{m}{Number of samples that are drawn.}
\item{n}{The size of each sample.}
\item{A}{The true slope value.}
\item{sig1}{Standard deviation of the measurement error for x.}
\item{sig2}{Standard deviation of the measurement error for y.}
\item{generator}{A function such that \code{generator(n)} gives the random vector of true values of X. }
}
\details{
    The vectors x=X+d and y=A*x+e are then the measured values where d,e are normally distributed with mean 0 and standard deviations sig1, sig2.  
    For these vectors the slope estimates by standard regression, orthogonal regression and case based latent regression are calculated. 
    }
\value{  A vector with the standard deviations of the slope estimates of standard regression, orthogonal regression and case based latent regression in this order.
}

\examples{
testRegs(100,50,0.6,0.2,0.4);
# this prints standard deviations of errors, e.g. 0.4014260 0.3478160 0.2742257 which shows that in this sampling CBregress performed best
    
# a more ambitous use is the following simulation study
print("Simulations with uniform distribution")
print(c("std reg.","Ortho.","Case based"))
for(A in c(0.5,1.0,1.5)) {
  res=c(0,0,0)
  print("")
  print(paste("Data simulation for a=",toString(A)))
  for(sigx in c(0.2,0.4)) for(sigy in c(0.2,0.4)) {
    print(paste("sig error x=",toString(sigx),"sig error y=",toString(sigy)))
    res=res+testRegs(10,100,A,sigx,sigy,generator=runif)
  } 
  print(c("Total errors",sqrt(res/4)))
}
    
}
