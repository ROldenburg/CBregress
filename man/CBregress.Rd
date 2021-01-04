\name{CBregress}
\alias{CBregress}
\title{Case based linear regression}
\description{
Fit a case based latent variable regression model.}
\usage{
CBregress(x,y)
}
\arguments{
\item{x,y}{Numeric vectors of equal length.}
}
\details{
    The function \code{CBregress} performs a case based univariate linear regression, similairly to \code{lm(y~x)}.
    In contrast to \code{lm}, the statistical model assumes and estimates a latent variable X such that x-X and y-(a*X+b) are as small errors as possible. 
}
\value{
  A list with the following fields.
  \item{a}{The estimate for the slope parameter}
  \item{b}{The estimate for the intercept parameter}
  \item{v0}{The estimate for the variance of X}
  \item{vd}{The estimate for the error variance of x-X_}
  \item{ve}{The estimate for the error variance of y-(a*X+b)}
  \item{X0}{The vector of estimatea for the true value X }
}
\references{
Reinhard Oldenburg (2020). Structural Equation Modeling
Comparing Two Approaches. The Mathematica Journal, URL
https://content.wolfram.com/uploads/sites/19/2020/12/Oldenburg.pdf .}

\examples{
vx=c(1,2,3,4)
vy=c(3,4,4,5)
res=CBregress(vx,vy)
res$a
res$b
}
