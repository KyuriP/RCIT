#' RCoT - tests whether x and y are conditionally independent given z. Calls RIT if z is empty.
#' @param x Random variable x.
#' @param y Random variable y.
#' @param z Random variable z.
#' @param approx Method for approximating the null distribution. Options include:
#' "lpd4," the Lindsay-Pilla-Basak method (default),
#' "gamma" for the Satterthwaite-Welch method,
#' "hbe" for the Hall-Buckley-Eagleson method,
#' "chi2" for a normalized chi-squared statistic,
#' "perm" for permutation testing (warning: this one is slow but recommended for small samples generally <500 )
#' @param num_f Number of features for conditioning set. Default is 25.
#' @param num_f2 Number of features for non-conditioning sets. Default is 5.
#' @param seed The seed for controlling random number generation. Use if you want to replicate results exactly. Default is NULL.
#' @return A list containing the p-value \code{p} and statistic \code{Sta}
#' @export
#' @examples
#' RCIT(rnorm(1000),rnorm(1000),rnorm(1000));
#'
#' x=rnorm(10000);
#' y=(x+rnorm(10000))^2;
#' z=rnorm(10000);
#' RCIT(x,y,z,seed=2);


RCoT <- function(x, y, z = NULL, suffStat, approx = "lpd4", num_f = 100, num_f2 = 5, seed = NULL) {
  
  # Extract data columns for x, y, and z from suffStat$data based on indices
  x_data <- suffStat$data[, x]
  y_data <- suffStat$data[, y]
  
  if (!is.null(z)) {
    z_data <- suffStat$data[, z]
  } else {
    z_data <- NULL
  }
  
  # If z is empty, call RIT for unconditional independence
  if (length(z_data) == 0) {
    out <- RIT(x_data, y_data, approx = approx, seed = seed)
    return(out)
  } else {
    # # Handle conditional independence with RCoT
    # x_data <- matrix2(x_data)
    # y_data <- matrix2(y_data)
    # z_data <- matrix2(z_data)
    
    # Filter out columns of z with zero standard deviation
    z_data <- z_data[, apply(z_data, 2, sd) > 0]
    # z_data <- matrix2(z_data)
    d <- ncol(z_data)
    
    if (length(z_data) == 0) {
      out <- RIT(x_data, y_data, approx = approx, seed = seed)
      return(out)
    } else if (sd(x_data) == 0 | sd(y_data) == 0) {
      # return(list(p = 1, Sta = 0))
      return(p = 1) # just pval
    }
    
    r <- nrow(x_data)
    if (r > 500) {
      r1 <- 500
    } else {
      r1 <- r
    }
    
    x_data <- normalize(x_data)
    y_data <- normalize(y_data)
    z_data <- normalize(z_data)
    
    # Generate random Fourier features for x, y, and z
    four_z <- random_fourier_features(z_data[, 1:d], num_f = num_f, sigma = median(c(t(dist(z_data[1:r1, ])))), seed = seed)
    four_x <- random_fourier_features(x_data, num_f = num_f2, sigma = median(c(t(dist(x_data[1:r1, ])))), seed = seed)
    four_y <- random_fourier_features(y_data, num_f = num_f2, sigma = median(c(t(dist(y_data[1:r1, ])))), seed = seed)
    
    # Normalize the random Fourier features
    f_x <- normalize(four_x$feat)
    f_y <- normalize(four_y$feat)
    f_z <- normalize(four_z$feat)
    
    # Compute covariance matrices
    Cxy <- cov(f_x, f_y)
    Czz <- cov(f_z)
    
    # Compute inverse of Czz using Cholesky decomposition
    i_Czz <- chol2inv(chol(Czz + diag(num_f) * 1e-10))
    Cxz <- cov(f_x, f_z)
    Czy <- cov(f_z, f_y)
    
    # Compute residuals
    z_i_Czz <- f_z %*% i_Czz
    e_x_z <- z_i_Czz %*% t(Cxz)
    e_y_z <- z_i_Czz %*% Czy
    
    res_x <- f_x - e_x_z
    res_y <- f_y - e_y_z
    
    # Approximate null distributions
    if (num_f2 == 1) {
      approx <- "hbe"
    }
    
    if (approx == "perm") {
      Cxy_z <- cov(res_x, res_y)
      Sta <- r * sum(Cxy_z^2)
      
      nperm <- 1000
      Stas <- c()
      for (ps in 1:nperm) {
        perm <- sample(1:r, r)
        Sta_p <- Sta_perm(res_x[perm, ], res_y, r)
        Stas <- c(Stas, Sta_p)
      }
      
      p <- 1 - (sum(Sta >= Stas) / length(Stas))
    } else {
      Cxy_z <- Cxy - Cxz %*% i_Czz %*% Czy
      Sta <- r * sum(Cxy_z^2)
      
      d <- expand.grid(1:ncol(f_x), 1:ncol(f_y))
      res <- res_x[, d[, 1]] * res_y[, d[, 2]]
      Cov <- 1/r * (t(res) %*% res)
      
      if (approx == "chi2") {
        i_Cov <- ginv(Cov)
        Sta <- r * (c(Cxy_z) %*% i_Cov %*% c(Cxy_z))
        p <- 1 - pchisq(Sta, length(c(Cxy_z)))
      } else {
        eig_d <- eigen(Cov, symmetric = TRUE)
        eig_d$values <- eig_d$values[eig_d$values > 0]
        
        if (approx == "gamma") {
          p <- 1 - sw(eig_d$values, Sta)
        } else if (approx == "hbe") {
          p <- 1 - hbe(eig_d$values, Sta)
        } else if (approx == "lpd4") {
          eig_d_values <- eig_d$values
          p <- try(1 - lpb4(eig_d_values, Sta), silent = TRUE)
          if (!is.numeric(p) | is.nan(p)) {
            p <- 1 - hbe(eig_d$values, Sta)
          }
        }
      }
    }
    
    if (p < 0) p <- 0
    
    # out=list(p=p,Sta=Sta);
    out = p # extract pval only
    return(out)
    }
}
 