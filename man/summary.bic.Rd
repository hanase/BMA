\name{summary.bic}
\alias{summary}
\alias{summary.bicreg}
\alias{summary.bic.glm}
\alias{summary.bic.surv}
\alias{summary.glib}
\alias{summary.mc3}
\alias{print}
\alias{print.bicreg}
\alias{print.bic.glm}
\alias{print.bic.surv}
\alias{print.glib}
\alias{print.mc3}

\title{Summaries of Bayesian model averaging objects}
\description{
     \code{summary} and \code{print} methods for Bayesian model averaging objects.
}
\usage{

\method{summary}{bicreg}(object, n.models = 5, 
         digits = max(3, getOption("digits") - 3), 
         conditional = FALSE, display.dropped = FALSE, ...)

\method{summary}{bic.glm}(object, n.models = 5, 
         digits = max(3, getOption("digits") - 3), 
         conditional = FALSE, display.dropped = FALSE, ...)

\method{summary}{bic.surv}(object, n.models = 5, 
         digits = max(3, getOption("digits") - 3), 
         conditional = FALSE, display.dropped = FALSE, ...)

\method{summary}{glib}(object, n.models = 5, 
         digits = max(3, getOption("digits") - 3), 
         conditional = FALSE, index.phi=1, ...) 

\method{summary}{mc3}(object, n.models = 5, 
         digits = max(3, getOption("digits") - 3), ...)

\method{print}{bicreg}(x, digits = max(3, getOption("digits") - 3), ...)

\method{print}{bic.glm}(x, digits = max(3, getOption("digits") - 3), ...)

\method{print}{bic.surv}(x, digits = max(3, getOption("digits") - 3), ...)

\method{print}{mc3}(x, digits = max(3, getOption("digits") - 3), 
                    n.models = nrow(x$variables), ...)

}
\arguments{
  \item{object}{ object of type 'bicreg', 'bic.glm', 'bic.surv', 'glib' or 'mc3'}
  \item{x}{ object of type 'bicreg', 'bic.glm', 'bic.surv', 'glib' or 'mc3'}
  \item{n.models}{ optional number specifying the number of models to display in summary }
  \item{digits}{ optional number specifying the number of digits to display }
  \item{conditional}{ optional logical value specifying whether to display conditional expectation and standard deviation}
  \item{display.dropped}{ optional logical value specifying whether to display the names of any variables dropped before model averaging takes place }
  \item{index.phi}{ optional number specifying which value of phi to use if multiple values of phi were run. Applies to \code{glib} objects only}
  \item{...}{other parameters to be passed to \code{print.default}}
}
\details{
The print methods display a view similar to \code{print.lm} or \code{print.glm}. 
The summary methods display a view specific to model averaging. 
}
\author{Ian Painter \email{ian.painter@AT@gmail.com}}
\note{The summary function does not create a summary object (unlike \code{summary.lm} or \code{summary.glm}), instead it directly prints the summary. Note that no calculations are done to create the summary.}
\examples{
# logistic regression
library("MASS")
data(birthwt)
y<- birthwt$lo
x<- data.frame(birthwt[,-1])
x$race<- as.factor(x$race)
x$ht<- (x$ht>=1)+0
x<- x[,-9]
x$smoke <- as.factor(x$smoke)
x$ptl<- as.factor(x$ptl)
x$ht <- as.factor(x$ht)
x$ui <- as.factor(x$ui)

glm.out1<- bic.glm(x, y, OR = 20, glm.family="binomial", 
                   factor.type=TRUE)
summary(glm.out1, conditional = TRUE)
}
\keyword{print}
