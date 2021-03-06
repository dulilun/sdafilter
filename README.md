
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

## Example 1: multiple testing: one-sample t-test

This is a simple example of multiple testing of means under dependence.

``` r
library(sdafilter)
n = 50
p = 100

rho = 0.8
Sig = matrix(rho, p, p)
diag(Sig) = 1

dat <- MASS::mvrnorm(n, rep(0, p), Sig)
mu = rep(0, p)
mu[1:as.integer(0.1*p)]=0.5
dat = dat+rep(1, n)%*%t(mu)

alpha = 0.2
#kwd = {'lasso', 'de_lasso', 'innovate', 'pfa'}
out = SDA_M(dat, alpha, solve(Sig), kwd='innovate')
print(out)
#>  [1]  1  2  3  4  5  6  7  8  9 10 15 54
```

## Example 2: multiple testing: two-sample t-test

This is a demonstration on how to use our SDA method in two sample case.

``` r
p = 100
n = 30

dat_I = matrix(rnorm(n*p),nrow = n)
mu = rep(0, p)
mu[1:10] = 1.5
dat_I = dat_I = rep(1, n)%*%t(mu)

dat_II = matrix(rnorm(n*p), nrow = n)

Sigma_I = diag(p)
Sigma_II = diag(p)

out = SDA_2S(dat_I, dat_II, alpha=0.05, Sigma_I, Sigma_II)
print(out)
#>  [1]  1  2  3  4  5  6  7  8  9 10
```
