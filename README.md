
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sdafilter

<!-- badges: start -->

<!-- badges: end -->

The goal of sdafilter is to provide an R package for our new paper on
arXiv; see the link: <https://arxiv.org/abs/2002.11992>

## Installation

You can install the released version of sdafilter from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sdafilter")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dulilun/sdafilter")
```

## Example

This is a basic example which shows you how to solve a simple problem:

``` r
library(sdafilter)
n = 50
p = 100
dat = matrix(rnorm(n*p), nrow=n)
mu = rep(0, p)
mu[1:as.integer(0.1*p)]=0.3
dat = dat+rep(1, n)%*%t(mu)
alpha = 0.2
out = SDA_M(dat, alpha, diag(p))
#> Warning in min(t[which(Ta <= alpha)]): no non-missing arguments to min;
#> returning Inf

#> Warning in min(t[which(Ta <= alpha)]): no non-missing arguments to min;
#> returning Inf
print(out)
#> [1] 5 7
```
