\name{MC3.REG.choose}
\alias{MC3.REG.choose}
\title{ Helper function to MC3.REG }
\description{
  Helper function to MC3.REG that chooses the proposal model for a Metropolis-Hastings step.
}
\usage{
MC3.REG.choose(M0.var, M0.out)
}
\arguments{
\item{M0.var}{a logical vector specifying the variables in the current model.}
\item{M0.out}{a logical vector specifying the outliers in the current model.}
}
\value{
  A list representing the proposal model, with components
  \item{var}{a logical vector specifying the variables in the proposal model.}
  \item{out}{a logical vector specifying the outliers in the proposal model.}
}
\references{    


Bayesian Model Averaging for Linear Regression Models 
Adrian E. Raftery, David Madigan, and Jennifer A. Hoeting (1997). 
Journal of the American Statistical Association, 92, 179-191. 

A Method for Simultaneous Variable and Transformation Selection in Linear Regression 
Jennifer Hoeting, Adrian E. Raftery and David Madigan (2002). 
Journal of Computational and Graphical Statistics 11 (485-507) 

A Method for Simultaneous Variable Selection and Outlier Identification in Linear Regression 
Jennifer Hoeting, Adrian E. Raftery and David Madigan (1996). 
Computational Statistics and Data Analysis, 22, 251-270 

Earlier versions of these papers are available via the World Wide Web using the url:
\url{http://www.stat.colostate.edu/~jah/papers/}

}
\author{Jennifer Hoeting \email{jah@AT@stat.colostate.edu} with the assistance of Gary Gadbury. Translation from Splus to R by Ian Painter \email{ian.painter@AT@gmail.com}.}
\note{ The implementation here differs from the Splus implentation. The Splus implementation uses global variables to contain the state of the current model and the history of the Markov-Chain. This implentation passes the current state and history to the function and then returns the updated state. }
\seealso{ \code{\link{MC3.REG}}, \code{\link{For.MC3.REG}}, \code{\link{MC3.REG.logpost}} }
\keyword{regression}
