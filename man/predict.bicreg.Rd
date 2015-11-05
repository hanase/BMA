\name{predict.bicreg}
\alias{predict.bicreg}
\title{Predict function for Bayesian Model Averaging for linear models.}
\description{
Bayesian Model Averaging (BMA) accounts for the model uncertainty inherent in the 
variable selection problem by averaging over the best models in the model class 
according to approximate posterior model probability. This function predicts
the response resulting from a BMA linear model from given data.
}
\usage{
\method{predict}{bicreg}( object, newdata, quantiles, \dots)
}
\arguments{
  \item{object}{a fitted object inheriting from class \code{bicreg}.}
  \item{newdata}{a data frame containing observations on variables from 
                 which the predictor variables are to be selected or 
                 constructed from a formula.}
  \item{quantiles}{The quantiles for which a predictive estimate is
    desired. The default is \code{c(.1,.5,.9)}, corresponding to the
    median (.5), and the 10th and 90th precentiles.}
  \item{\dots}{ignored (for compatibility with generic function).}
}
\value{
  The predicted response values from the BMA model for each observation in newdata.
}

\seealso{ \code{\link{bicreg}}  }
\examples{

  library(MASS)

# Example 1

     data(UScrime)

     x <- UScrime[,-16]
     y <- log(UScrime[,16])
     x[,-2]<- log(x[,-2])

     crimeBMA <- bicreg(x, y, strict = FALSE, OR = 20)
     predict( crimeBMA, x)

# Example 2 (Venables and Ripley)

     npkBMA <- bicreg( x = npk[, c("block","N","K")], y=npk$yield)
     predict( npkBMA, newdata = npk)

# Example 2 (Venables and Ripley)

     gasPRbma <- bicreg( x = whiteside[,c("Insul", "Temp")], 
                         y = whiteside$Gas)
     predict( gasPRbma, newdata = whiteside)

}
\keyword{regression}
\keyword{models}
