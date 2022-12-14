% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irtpwr.R
\name{irtpwr}
\alias{irtpwr}
\title{Perform Power Analysis}
\usage{
irtpwr(
  hyp,
  stat = c("Wald", "LR", "Score", "Gradient"),
  method = "analytical",
  sampling.npers = 10^5,
  approx.npers = 10^5,
  SE.type = "Oakes",
  sampling.mat = "ApproxFisher",
  power = NULL,
  N = NULL,
  alpha = NULL
)
}
\arguments{
\item{hyp}{Hypothesis Object created by the setup.hypothesis function}

\item{stat}{character vector containing the statistics to be calculated. Options are 'Wald','LR','Score', and 'Gradient'. By default, all statistics are included}

\item{method}{character, indicating the method used. The options are 'analytical'(default) for the analytical power analysis method or 'sampling' for the sampling-based method. The sampling-based method is generally recommended for higher numbers of items.}

\item{sampling.npers}{integer, sample size for the sampling-based approach. An artificial data set of this size is generated to fit a model and later estimate the noncentrality parameter from.}

\item{approx.npers}{integer, sample size for approximating the Fisher expected information matrix in the sampling-based approach. An artificial data set is calculated of this size to calculate the Fisher expected information matrix from. In contrast to the data created with the sampling.npers sample size, this sample is not used to fit a model.}

\item{SE.type}{Method for calculation of the observed information matrix used for calculating the statistics in the sampling-based approach ('Oakes' by default). Another option is 'Fisher'.}

\item{sampling.mat}{Approach to calculate the information matrix used for calculating the statistics in the sampling-based approach. By default ('ApproxFisher'), an sampling-based approximation of the expected Fisher matrix is calculated using an observed information matrix of the type SE.type}

\item{power}{numeric, statistical power for which the necessary sample size is calculated}

\item{N}{integer, sample size for which the statistical power is calculated.}

\item{alpha}{numeric, alpha level}
}
\value{
function returns an object of class irtpwrresult
}
\description{
Perform analytical or sampling-based power analysis for the Wald, LR, score, or gradient statistic.
}
\examples{
\donttest{

library(mirt)
dat <- expand.table(LSAT7)
mirtfit <- mirt(dat,1,verbose = FALSE)
hyp <- setup.hypothesis(type = '1PLvs2PL', altpars = mirtfit)
res <- irtpwr(hyp=hyp,alpha=.05,power =.8)
summary(res)
}
}
