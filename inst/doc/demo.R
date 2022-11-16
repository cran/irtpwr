## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "70%",
  fig.dim = c(6, 4),
  # tidy = TRUE,
  # tidy.opts=list(arrow=TRUE,width.cutoff = 50),
  eval=T
)

## ----include=FALSE------------------------------------------------------------
set.seed(1)

## -----------------------------------------------------------------------------
library(irtpwr)
library(mirt)

## -----------------------------------------------------------------------------
dat <- expand.table(LSAT7)
mirtfit <- mirt(dat, 1, verbose = FALSE)

## -----------------------------------------------------------------------------
hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)

## -----------------------------------------------------------------------------
res <- irtpwr(hyp = hyp, power = 0.8, alpha = 0.05)
summary(res)

## -----------------------------------------------------------------------------
plot(res)

## -----------------------------------------------------------------------------
altpars <- list(a = rlnorm(5, sdlog = 0.4), d = rnorm(5))
hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = altpars)

## -----------------------------------------------------------------------------
group1 <- group2 <- list(a = rlnorm(5, sdlog = 0.2),
    d = rnorm(5))

group2$a[1] <- (group2$a[1])^2
group2$d[1] <- group2$d[1] + 0.5

altpars <- list(group1, group2)

hyp <- setup.hypothesis(type = "DIF2PL", altpars = altpars)


## -----------------------------------------------------------------------------
res <- irtpwr(hyp = hyp, N = 600, alpha = 0.05)
summary(res)

## -----------------------------------------------------------------------------
summary(res, power = 0.8)

## -----------------------------------------------------------------------------
summary(res, N = 700)

## -----------------------------------------------------------------------------
res <- irtpwr(hyp = hyp, N = 600, alpha = 0.05, method = "sampling")
summary(res)

## -----------------------------------------------------------------------------
calc.time(hyp, n.items = 7)

## -----------------------------------------------------------------------------
altpars <- list(
  a = c(1.1,seq(.5,1.4,length.out=4)),
  d = c(.1,seq(1.3,-1.3,length.out=4)),
  g = c(.22,rep(.2,4))
)

## -----------------------------------------------------------------------------
hyp <- setup.hypothesis(type = "3PL_basic", altpars = altpars)

## -----------------------------------------------------------------------------
res <- irtpwr(hyp = hyp, alpha = 0.05, power = 0.8,
    method = "sampling", SE.type = "Fisher")

## -----------------------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
dat <- expand.table(LSAT7)
dat <- dat[, c(2, 4, 1, 3, 5)]  # re-ordering items so that items 1 and 2 of the resulting data frame are compared
model_an <- "F1 = 1-5
      F2 = 1-4"
altpars <- mirt(dat, model = mirt.model(model_an),
    verbose = FALSE)

## -----------------------------------------------------------------------------
hyp <- setup.hypothesis(type = "multi_basic", altpars = altpars)

## -----------------------------------------------------------------------------
res <- irtpwr(hyp = hyp, alpha = 0.05, power = 0.8,
    method = "sampling", SE.type = "Fisher")

## -----------------------------------------------------------------------------
summary(res)

