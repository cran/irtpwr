## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
    out.width = "70%",
  # tidy = TRUE,
  # tidy.opts=list(arrow=TRUE,width.cutoff = 50),
  eval=T
)

## ----message=FALSE, echo = FALSE----------------------------------------------
library(irtpwr)

## -----------------------------------------------------------------------------
res <- function(altpars, nullpars = NULL) {

    n.items <- length(altpars[[1]]) # we can read off the number of items from the altpars object
    
    # the A matrix represents the calculations that need to be performed on the item parameters according to the hypothesis. in This case, only the difficulty of the first item needs to be extracted. 
    Amat <- c(0, 1, rep(0, (n.items - 1) * 2)) |>
            (function(x) matrix(x, ncol = n.items *
                2, byrow = TRUE))()
    # the c vector is the value that the item parameters are compared to after transformation by the A matrix. In this case, the difficulty parameter is only compared against 0. 
    cvec = 0
    # By specifying a mirt.model, we instruct mirt on how to fit the restricted model. In this case, it is a model where the first difficulty parameter is kept at 0. 
    model = mirt::mirt.model(paste("F = 1-",
            n.items, "
                           FIXED = (1, d)
                           START = (1,d,0)"))

    re <- list(n.items = n.items, itemtype = "2PL",
        Amat = Amat, cvec = cvec, model = model)
    return(re)
}

## -----------------------------------------------------------------------------
unres <- function(altpars) {

   # We first transform the parameters altpars from a list to a vector using the longpars argument. It results from a concatenation of the discrimination and difficulty parameters.
    longpars = pars.long(pars = altpars, itemtype = "2PL")
    # For the unrestricted model, we fit a simple 2PL model. This is specified in mirt using a 1 and the "2PL" itemtype. 
    model = 1
    itemtype = "2PL"
    re <- list(parsets = altpars, model = model, itemtype = itemtype, longpars = longpars)

    return(re)
}

## -----------------------------------------------------------------------------
maximizeL <- function(hyp) {
    # Hypothesis-specific algorithm to find the
    # maximum likelihood restricted parameter set
  
  # in this case, the a parameter of the first item is searched for.


    maxlpreload <- function(pars) {
        # returns the density for each response
        # pattern under the model parameters pars

      # setting up all response patterns 
        patterns <- as.matrix(expand.grid(lapply(seq_len(length(pars$a)),
            function(x) c(0, 1))))


        pre <- c()
        #calculating the density
        for (i in seq_len(nrow(patterns))) {
            pre[i] <- funs$g(patterns[i, ], pars)
        }

        return(pre)
    }


    maxl <- function(x, pars, pre) {
        # calculates the likelihood of parameters
        # x given model 'pars'
      
      # setting up all response patterns 
        patterns <- as.matrix(expand.grid(lapply(seq_len(length(pars$a)),
            function(x) c(0, 1))))

        # collecting all parameters in a list, inserting x as the first a parameter
        x <- list(a = c(x, pars$a[2:length(pars$a)]),
            d = c(0, pars$d[2:length(pars$d)]))

        res <- c()
        # calculating the likelihoods for each pattern under this parameter set
        for (i in seq_len(nrow(patterns))) {
            px <- pre[i]
            qx <- funs$g(patterns[i, ], x)
            res[i] <- {
                px * log(qx)
            }
        }
        re <- -sum(res)
    }
    resmod <- hyp$resmod
    unresmod <- hyp$unresmod

    pars <- unresmod$parsets
    # loading the model specific density functions 
    funs <- load.functions(unresmod$itemtype)

    # setting some starting value for the optimization
    startval <- pars$a[1]

    # calculating the densities as definied above
    maxlpre <- maxlpreload(pars)

    # finding the a parameter with the highest likelihood
    optpar <- stats::optim(startval, function(x) {
        maxl(x, pars, maxlpre)
    }, method = "BFGS")
    re <- pars
    # saving the resulting item parameters
    re$a <- c(optpar$par[1], pars$a[2:length(pars$a)])
    re$d <- c(0, pars$d[2:length(pars$d)])

    return(re)
}


## -----------------------------------------------------------------------------
h_2PL_basic <- list(res = res, unres = unres, maximizeL = maximizeL)

## -----------------------------------------------------------------------------
altpars <- list(a = rlnorm(5, sdlog = 0.4), d = rnorm(5))

altpars$d[1] <- 0.2

hyp <- setup.hypothesis(type = h_2PL_basic, altpars = altpars)

res <- irtpwr(hyp = hyp, alpha = 0.05, power = 0.8)
summary(res)

