## ----setup.hypothesis usage, results = "hide", eval = FALSE-------------------
#  setup.hypothesis(type, altpars = NULL, nullpars = NULL)

## ----install packages, results = "hide", eval = FALSE-------------------------
#  install.packages("irtpwr")

## ----load packages, results = "hide", eval = TRUE, message = FALSE------------
library("irtpwr")

## ----set.seed Rasch vs 2PL a priori, echo = FALSE-----------------------------
set.seed(7443)

## ----Rasch vs 2PL setup.hypothesis a priori, results = "hide", eval = TRUE----
# parameters under the alternative hypothesis
altpars <- list(a = rlnorm(5, sdlog = 0.4),    
                d = rnorm(5))                  
 
hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = altpars)

## ----Rasch vs 2PL irtpwr a priori, results = "hide", eval = TRUE--------------
res <- irtpwr(hyp = hyp, power = 0.8, alpha = 0.05)
summary(res)

#> Sample sizes for power = 0.8 (alpha = 0.05): 
#>
#> Statistic   N
#>      Wald 503
#>        LR 504
#>     Score 520
#>  Gradient 505
#>
#> Method: Analytical


## ----Rasch vs 2PL irtpwr a priori plot, echo = TRUE---------------------------
plot(res) 

## ----Rasch vs 2PL irtpwr a priori pwr 0.9, results = "hide", eval = TRUE------
summary(res, power = 0.9)

#> Sample sizes for power = 0.8 (alpha = 0.05): 
#>
#> Statistic   N
#>      Wald 649
#>        LR 650
#>     Score 671
#>  Gradient 652
#>
#> Method: Analytical


## ----install packages mirt, results = "hide", eval = FALSE--------------------
#  install.packages("mirt")

## ----load packages mirt, results = "hide", eval = TRUE, message = FALSE-------
library("mirt")

## ----set.seed Rasch vs 2PL a posteriori, echo = FALSE-------------------------
set.seed(1258)

## ----Rasch vs 2PL data a posteriori, results = "hide", eval = TRUE------------
dat <- expand.table(LSAT7)

## ----Rasch vs 2PL mirt a posteriori, results = "hide", eval = TRUE------------
mirtfit <- mirt(dat, 1, verbose = FALSE)

## ----Rasch vs 2PL setup.hypothesis a posteriori, results = "hide", eval = TRUE----
hyp <- setup.hypothesis(type = "1PLvs2PL", altpars = mirtfit)

## ----Rasch vs 2PL irtpwr a posteriori, results = "hide", eval = TRUE----------
res <- irtpwr(hyp = hyp, N = 1000, alpha = 0.05)
summary(res)

#> Power for N = 1000 (alpha = 0.05):  
#>
#> Statistic  Power
#>      Wald 0.7411
#>        LR 0.8513
#>     Score 0.8308
#>  Gradient 0.8987
#>
#> Method: Analytical


## ----Rasch vs 2PL irtpwr a posteriori plot, echo = TRUE-----------------------
plot(res) 

## ----Rasch vs 2PL irtpwr sampling a posteriori, results = "hide", eval = TRUE----
res <- irtpwr(hyp = hyp, N = 600, alpha = 0.05, method = "sampling")
summary(res)

#> Power for N = 1000 (alpha = 0.05):  
#>
#> Statistic  Power
#>      Wald 0.4946
#>        LR 0.6073
#>     Score 0.5957
#>  Gradient 0.6578
#>
#> Method: Sampling-Based


## ----calc.time, results = "hide", eval = TRUE---------------------------------
calc.time(hyp, n.items = 10)

#> [1] 103.488

## ----set.seed, echo = FALSE---------------------------------------------------
set.seed(8989)

## ----DIF setup.hypothesis a priori, results = "hide", eval = TRUE-------------
native <- nonnative <- list(a = rlnorm(5, sdlog = 0.2), 
                            d = rnorm(5))               

nonnative$a[1] <- (nonnative$a[1])^2   # DIF effect in the a parameter
nonnative$d[1] <- nonnative$d[1] - 0.5 # DIF effect in the d parameter

altpars <- list(native, nonnative) # parameters under the alternative hypothesis

hyp <- setup.hypothesis(type = "DIF2PL", altpars = altpars)

## ----DIF irtpwr a priori, results = "hide", eval = TRUE-----------------------
res <- irtpwr(hyp = hyp, power = 0.8, alpha = 0.05)
summary(res)

#> Sample sizes for power = 0.8 (alpha = 0.05): 
#>
#> Statistic   N
#>      Wald 814
#>        LR 794
#>     Score 796
#>  Gradient 790
#>
#> Method: Analytical


## ----irtpwr a priori plot, echo = TRUE----------------------------------------
plot(res) 

## ----generate difdata, echo = FALSE-------------------------------------------
sim_2pl <- function(n_Per = 1000, n_Item = 5, covariate = "binary", a_delta = 1, d_delta = 0, t_delta = 0, numsplit = 50) {
  a <- runif(n_Item, 1, 2.25)
  d <- rnorm(n_Item, 0, 1)
  a_diff <- a
  d_diff <- d
  # Kein DIF auf Item 1 | a_delta = 1, d_delta = 0
  a_diff[1] <- a_diff[1] * a_delta
  d_diff[1] <- d_diff[1] - d_delta
  if (covariate == "binary") {
    language <- sample(c("native", "nonnative"), n_Per, replace = TRUE, prob = c(.5, .5))
    language <- sort(language)
    t_ref <- rnorm(sum(language == "native"), 0 - t_delta / 2, 1)
    t_dif <- rnorm(sum(language == "nonnative"), 0 + t_delta / 2, 1)
    dat1 <- simdata(a, d, N = sum(language == "native"), itemtype = "2PL", Theta = t_ref)
    dat2 <- simdata(a_diff, d_diff, N = sum(language == "nonnative"), itemtype = "2PL", Theta = t_dif)
    dat3 <- rbind(dat1, dat2)
    colnames(dat3) <- NULL
    dat <- data.frame(matrix(nrow = n_Per, ncol = 0))
    dat$item <- as.matrix(dat3)
    dat$language <- as.factor(language)
    return(dat)
  } else if (covariate == "numeric") {
    age <- round(runif(n_Per, 0, 100))
    while (median(age) != 50) {
      age <- round(runif(n_Per, 0, 100))
    }
    age <- sort(age)
    obs <- which(age == numsplit)
    t_ref <- rnorm(obs[length(obs)], 0 - t_delta / 2, 1)
    t_dif <- rnorm(n_Per - obs[length(obs)], 0 + t_delta / 2, 1)
    dat1 <- simdata(a, d, N = obs[length(obs)], itemtype = "2PL", Theta = t_ref)
    dat2 <- simdata(a_diff, d_diff, N = n_Per - obs[length(obs)], itemtype = "2PL", Theta = t_dif)
    dat3 <- rbind(dat1, dat2)
    colnames(dat3) <- NULL
    dat <- data.frame(matrix(nrow = n_Per, ncol = 0))
    dat$Item <- as.matrix(dat3)
    dat$age <- age
    return(dat)
  } else {
    stop("Covariate is neither binary nor numeric")
  }
}

# generate data
set.seed(5666)
dat <- sim_2pl(732, 5, "binary", a_delta = 2, d_delta = 0.5, t_delta = 0)

## ----showing data, results = "hide"-------------------------------------------
head(dat)

#>   item.1 item.2 item.3 item.4 item.5  language
#> 1      0      0      0      0      0 nonnative
#> 2      0      0      0      0      1    native
#> 3      1      1      0      1      1    native
#> 4      0      0      1      0      0 nonnative
#> 5      1      0      1      0      1    native
#> 6      0      0      0      1      1 nonnative

## ----DIF mirt a posteriori, results = "hide", eval = TRUE---------------------
newmodel <- '
    F = 1-5
    CONSTRAINB = (2-5, a1), (2-5, d)'

difmirt <- multipleGroup(as.data.frame(dat[,"item"]), newmodel, group = dat$language)

## ----DIF coef, results = "hide", eval = TRUE----------------------------------
pars <- coef(difmirt, 
     IRTpars = FALSE,  
     simplify = TRUE)
pars

#> $native
#> $items
#>       a1      d
#> V1 1.558  0.594
#> V2 2.241 -0.390
#> V3 2.204 -0.328
#> V4 1.908 -1.043
#> V5 1.003  1.040
#> 
#> $nonnative
#> $items
#>       a1      d
#> V1 1.677  0.012
#> V2 2.241 -0.390
#> V3 2.204 -0.328
#> V4 1.908 -1.043
#> V5 1.003  1.040

## ----DIF setup.hypothesis a posteriori, results = "hide", eval = TRUE---------
# parameters of native English speaking test taker
pars_native <- pars[["native"]][["items"]]
native <- list(a = pars_native[,"a1"], 
               d = pars_native[,"d"])

# parameters of non-native English speaking test taker
pars_nonnative <- pars[["nonnative"]][["items"]]
nonnative <- list(a = pars_nonnative[,"a1"], 
                  d = pars_nonnative[,"d"])

altpars <- list(native, nonnative) # parameters under the alternative hypothesis

hyp <- setup.hypothesis(type = "DIF2PL", altpars = altpars)

## ----DIF irtpwr a posteriori, results = "hide", eval = TRUE-------------------
res <- irtpwr(hyp = hyp, N = 732, alpha = 0.05)
summary(res)

#>  Power for N = 732 (alpha = 0.05): 
#> 
#>  Statistic  Power
#>       Wald 0.7779
#>         LR 0.7917
#>      Score 0.7900
#>   Gradient 0.7949
#> 
#> Method: Analytical

## ----irtpwr a posteriori plot, echo = TRUE------------------------------------
plot(res) 

## ----irtpwr pwr 0.9, results = "hide", eval = TRUE----------------------------
summary(res, power = 0.9) 

#>  Sample sizes for power = 0.9 (alpha = 0.05): 
#> 
#>  Statistic    N
#>       Wald 1013
#>         LR  981
#>      Score  985
#>   Gradient  974
#> 
#> Method: Analytical

## ----irtpwr N 800, results = "hide", eval = TRUE------------------------------
summary(res, N = 800) 

#>  Power for N = 800 (alpha = 0.05): 
#> 
#>  Statistic  Power
#>       Wald 0.8154
#>         LR 0.8283
#>      Score 0.8267
#>   Gradient 0.8313
#> 
#> Method: Analytical

