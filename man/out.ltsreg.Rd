\name{out.ltsreg}
\alias{out.ltsreg}
\title{out.ltsreg }
\description{
 Function to identify potential outliers
}
\usage{
out.ltsreg(x, y, delta)
}
\arguments{
  \item{x}{the design matrix}
  \item{y}{observations}
  \item{delta}{ the threshold set by the user.  Standardized residuals from least trimmed squares regression that are larger than delta are identified as potential outliers}
}
\value{
  A 0/1 vector indicating whether each observation is a potential outlier. The function was designed for use with the variable and outlier selection function \code{MC3.REG}
}
\author{Jennifer A. Hoeting}
\seealso{  \code{\link{MC3.REG}}}
\keyword{regression}
