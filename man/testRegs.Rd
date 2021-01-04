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
\item{sig1}{Standard deviation of the measurement error for _x_.}
\item{sig2}{Standard deviation of the measurement error for _y_.}
\item{generator}{A function such that \code{generator(n)} gives the random vector of true values of _X_. }
}
\details{
    The vectors _x=X+d_ and _y=A*x+e_ are then the measured values where _d,e_ are normally distributed with mean 0 and standard deviations _sig1, sig2_.  
    For these vectors the slope estimates by standard regression, orthogonal regression and case based latent regression are calculated. 
    }
\value{  A vector with the standard deviations of the slope estimates of standard regression, orthogonal regression and case based latent regression in this order.
}

\examples{
testRegs(100,50,0.6,0.2,0.4)
# this prints standard deviations of errors, e.g. 0.4014260 0.3478160 0.2742257 which shows that in this sampling CBregress performed best
}