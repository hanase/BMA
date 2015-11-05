\name{imageplot.bma}
\alias{imageplot.bma}
\title{Images of models used in Bayesian model averaging }
\description{
  Creates an image of the models selected using \code{\link{bicreg}}, \code{\link{bic.glm}} or \code{\link{bic.surv}}.
}
\usage{
imageplot.bma(bma.out, color = c("red", "blue", "#FFFFD5"), 
              order = c("input", "probne0", "mds"), ...) 
}
\arguments{
  \item{bma.out}{An object of type 'bicreg', 'bic.glm' or 'bic.surv'}
  \item{color}{A vector of colors of length 3, or a string with value "default" or "blackandwhite", representing the colors to use for the plot. 
  The first color is the color to use when the variable estimate is positive, the second color is the color to use when the variable estimate is negative, and the third color is the color to use when the variable is not included in the model.
  
  The value "default" is available for backward compatibility with the first version of \code{imageplot.bma}, and uses the same color for positive and negative estimates. 
  The value "blackandwhite" produces a black and white image.}
  \item{order}{The order in which to show the variables. The value "input" keeps the order as found in the object, the value "probne0" orders the variables in terms of probability of inclusion, and the value "mds" orders the variables using (single) multidimensional scaling} 
  \item{...}{Other parameters to be passed to the \code{image} and \code{axis} functions.}
}
\details{ 
Creates an image of the models selected using \code{\link{bicreg}}, \code{\link{bic.glm}} or \code{\link{bic.surv}}. The image displays inclusion and exclusion of variables within models using separate colors. By default the color for inclusion depends on whether the variable estimate for each model is positive or negative.

If the \code{factor.type == TRUE} option is set in the bma object being displayed, then \code{imageplot.bma} displays only inclusion and exclusion of models, with the color not linked to variable estimates.

The option \code{color = "mds"} is useful for observing variables with linked behavior, it attemps to order the variables in such a way as to keep variabiles with linked behavior (for example, one variabile is only included in a model when another variabile is not included in the model) close together.
This option uses multidimensional scaling on one dimension using Kendall's tau statistic calculated on two-by-two tables of pairwise comparisons of variable inclusion/exclusion from the selected models.
  }

\author{Adrian E. Raftery \email{raftery@AT@stat.washington.edu} and Hana Sevcikova}
\references{Clyde, M. (1999) Bayesian Model Averaging and Model Search Strategies (with discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid, J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages 157-185.
}
\seealso{ \code{\link{bicreg}}, \code{\link{bic.glm}}, \code{\link{bic.surv}} }
\examples{
# logistic regression using bic.glm
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

glm.out1<- bic.glm(x, y, strict = TRUE, OR = 20, glm.family="binomial")
imageplot.bma(glm.out1)

\dontrun{
# logistic regression using glib
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

glib.birthwt<- glib(x,y, error="binomial", link = "logit")
glm.birthwt<- as.bic.glm(glib.birthwt)
imageplot.bma(glm.birthwt, order = "mds")
}
}
\keyword{regression}
\keyword{models}
