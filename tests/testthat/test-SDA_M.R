test_that("SDA_M returns indices of rejected hypotheses", {

  n = 50
  p = 100
  dat = matrix(rnorm(n*p), nrow=n)
  mu = rep(0, p)
  mu[1:as.integer(0.1*p)]=0.5
  dat = dat+rep(1, n)%*%t(mu)
  alpha = 0.2

  out = SDA_M(dat, alpha, nonsparse = TRUE)
  if(is.character(out)){
    expect_match(out, 'no rejection')
  }else{
    a_in_b = all(out %in% 1:p)
    expect_equal(a_in_b, TRUE)
  }

})


test_that("test the kwd option in SDA_M", {

  n = 50
  p = 100

  # compound symmetry covariance
  rho = 0.8
  Sig = matrix(rho, p, p)
  diag(Sig) = 1

  dat <- MASS::mvrnorm(n, rep(0, p), Sig)
  mu = rep(0, p)
  mu[1:as.integer(0.1*p)]=0.5
  dat = dat+rep(1, n)%*%t(mu)
  alpha = 0.2
  out_1 = SDA_M(dat, alpha, nonsparse = TRUE, kwd='pfa')
  if(is.character(out_1)){
    expect_match(out_1, 'no rejection')
  }else{
    a_in_b = all(out_1 %in% 1:p)
    expect_equal(a_in_b, TRUE)
  }

})


test_that("test the kwd option in SDA_M again", {

  n = 50
  p = 100

  # AR(1) matrix
  ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                      (1:n - 1))
    rho^exponent
  }
  rho = 0.8
  Sig = ar1_cor(p, rho)

  dat <- MASS::mvrnorm(n, rep(0, p), Sig)
  mu = rep(0, p)
  mu[1:as.integer(0.1*p)]=0.5
  dat = dat+rep(1, n)%*%t(mu)
  alpha = 0.2
  out_1 = SDA_M(dat, alpha, solve(Sig), kwd='innovate')
  if(is.character(out_1)){
    expect_match(out_1, 'no rejection')
  }else{
    a_in_b = all(out_1 %in% 1:p)
    expect_equal(a_in_b, TRUE)
  }

})

