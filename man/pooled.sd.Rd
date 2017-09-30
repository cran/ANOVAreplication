\name{pooled.sd}
\alias{pooled.sd}
\title{Pooled standard deviation calculator}
\description{
Calculates the pooled standard deviation.}
\usage{
pooled.sd(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A dataframe with two variables: y for the dependent variable and g for the grouping variable. }
}
\value{
Returns the pooled standard deviation.
}
\author{
M.A.J. Zondervan-Zwijnenburg}
\examples{
data <- data.frame(y=ChickWeight$weight,g=ChickWeight$Diet)
pooled.sd(data)
}
