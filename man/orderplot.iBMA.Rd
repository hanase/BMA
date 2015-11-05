\name{orderplot}
\alias{orderplot}
\alias{orderplot.iBMA.glm}
\alias{orderplot.iBMA.bicreg}
\alias{orderplot.iBMA.surv}
\alias{orderplot.iBMA.intermediate.glm}
\alias{orderplot.iBMA.intermediate.bicreg}
\alias{orderplot.iBMA.intermediate.surv}

\title{ Orderplot of iBMA objects }
\description{
 This function displays a plot showing the selection and rejection of variables being considered in an iterated Bayesian model averaging variable selection procedure.
}

\usage{
orderplot(x, ...)

}
\arguments{
  \item{x}{ an object of type iBMA.glm, iBMA.bicreg, iBMA.surv, iBMA.intermediate.glm, iBMA.intermediate.bicreg or iBMA.intermediate.surv.}
  \item{\dots}{ other parameters to be passed to plot.default }
}
\details{
  The x-axis represents iterations, the y-axis variables. 
  For each variable, a dot in the far left indicates that the variable has not yet been examined, 
  a black line indicates the variable has been examined and dropped, the start of the line represents when the variable was first examined, 
  the end represents when the variable was dropped. A blue line represents a variable that is still in the selected set of variables.
  If the iterations have completed then the blue lines end with blue dots, representing the final set of variables selected.
  
}

\author{Ian Painter \email{ian.painter@AT@gmail.com}}
\seealso{ \code{\link{summary.iBMA.glm}}, \code{\link{iBMA}} }
\examples{

\dontrun{
############ iBMA.glm
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

### add 41 columns of noise
noise<- matrix(rnorm(41*nrow(x)), ncol=41)
colnames(noise)<- paste('noise', 1:41, sep='')
x<- cbind(x, noise)

iBMA.glm.out<- iBMA.glm(x, y,  glm.family="binomial", factor.type=FALSE, 
                        verbose = TRUE, thresProbne0 = 5 )
orderplot(iBMA.glm.out)
}

}
\keyword{ hplot }
