\name{Fbar.ineq}
\alias{Fbar.ineq}
\title{F-bar for inequality constraints}
%- Also NEED an '\alias' for EACH other topic documented here.
\description{
The function calculates F-bar for inequality constrained hypotheses (Type B). See Silvapulle & Sen (2011) for background on the F-bar statistic. The code of Vanbrabant (2017) is the basis to this Fbar function.
}
\usage{
Fbar.ineq(data,Amat)
}
\arguments{
  \item{data}{
A dataframe with two variables: (1) a dependent variable, and (2) a grouping variable.}
\item{Amat}{
A p by q matrix, where p is the number of means in the ANOVA model, and q is the number of constraints to be imposed on the model. Each row represents one constraint where the parameter with the lower value according to the constraint receives the value -1, and the parameter with the higher value according to the constraint receives the value 1. Other parameters within the same row obtain the value 0.}
}
\value{
The value for the F-bar statistic.
}
\references{
Silvapulle, M. J., & Sen, P. K. (2011). Constrained statistical inference: Order, inequality, and shape constraints (Vol. 912). John Wiley & Sons. doi: 10.1002/9781118165614.ch1

Vanbrabant, L. (2017). restriktor: Restricted Statistical Estimation and Inference for Linear Models. R package version 0.1-55. https://CRAN.R-project.org/package=restriktor
}
\author{
M. A. J. Zondervan-Zwijnenburg
}
\seealso{
See also \code{\link{runShiny}}, \code{\link{Fbar.dif}}, and \code{\link{Fbar.exact}}.
}
\examples{
data <- data.frame(y=ChickWeight$weight,g=ChickWeight$Diet)
aggregate(data$y,by=list(data$g),mean)

#make Amat with constraints: 1<4,2<4,3<4 (last constraint is not true)
Amat <- (rbind(c(-1,0,0,1),c(0,-1,0,1),c(0,0,-1,1)))

Fbar.ineq(data,Amat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

