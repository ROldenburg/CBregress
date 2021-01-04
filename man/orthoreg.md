\name{orthoreg}
\title{Linear univariate orthogonal regression}
\description{
Fit a line to two-dimensional data points by minimizing orthogonal ditances between line and points.}
\usage{
orthoreg(x,y)
}
\arguments{
\item{x,y}{Numeric vectors of equal length.}
}
\details{
    The function \code{orthoreg} performs an orthogonal linear regression. The model equation is _y=a*x+b_ 
}
\value{
  A vector with the estimates for _a_ and _b_. 
  }

\examples{
vx=c(1,2,3,4)
vy=c(3,4,4,5)
orthoreg(vx,vy)
}
