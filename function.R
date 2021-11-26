#-------------------------------------------------------------------------------
#                              Function: analysis
#-------------------------------------------------------------------------------
# Purpose: A function for penalized estimation of the single-index
#          varying-coefficient model. This function automatically generates grid
#          points for the tuning parameters and runs the entire solution path.
# Arguments:
# Y - Response variable. A Surv object for 'cox'.
# X - Matrix of covariates, of dimension n * p; each row is an observation
#     vector.
# U - Matrix of effect modifiers, of dimension n * q; each row is an observation
#     vector.
# Z - Matrix of unpenalized variables, of dimension n * r; each row is an
#     observation vector.
# nknot - Number of knots that define the spline that is the number of internal
#         breakpoints plus 2.
# degree - Non-negative integer degree of the piecewise polynomial.
# family - Response type of 'gaussian' or 'cox'.
# nlambda - The number of lambda values in each dimension.
# w - A vector of multiplicative factors by which each group's penalty is to be
#     multiplied.
# beta0 - A vector of initial estimate of the singe-index parameters.
# coef0 - A vector of initial estimate of the regression coefficients.
# lambda1.min.ratio - smallest value for lambda1, as a fraction of lambda.max.
# lambda2.min.ratio - smallest value for lambda2, as a fraction of lambda.max.
# maxit - Maximum number of replications. Default is 100.

analysis <- function(Y, X, U, Z, nknot, degree, family, nlambda = 20, w, beta0, 
                     coef0, lambda1.min.ratio, lambda2.min.ratio, maxit = 100) {
    
    n <- nrow(X)
    p <- ncol(X)
    q <- ncol(U)
    r <- ncol(Z)
    d <- nknot - 2 + degree
    
    if (missing(lambda1.min.ratio) & missing(lambda2.min.ratio)) {
        if (p == 20) {
            lambda1.min.ratio <- 0.001
            lambda2.min.ratio <- 0.001
        } else {
            lambda1.min.ratio <- 0.05
            lambda2.min.ratio <- 0.05
        }
    }
    
    # Generate knots based on U
    knots <- qnorm(seq(0, 1, l = nknot)[-c(1, nknot)])
    Boundary.knots <- c(- max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))), 
                        max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))))
    
    # Generate two dimensional grid points of lambda1 and lambda2
    lambda1 <- setup.grlambda(cbind(X, Z), Y, family = family, 
                              group = c(1:p, rep(0, r)), nlambda = nlambda, 
                              lambda.min.ratio = lambda1.min.ratio, 
                              penalty.factor = if (ncol(Z) == 0) w else c(w, 0))
    if (family == "cox") {
        fit.temp <- grpsurv(cbind(X, Z), Y, lambda = lambda1, 
                            group = c(1:p, rep(0, r)), max.iter = 1e+05, 
                            group.multiplier = w)
    }
    if (family == "gaussian") {
        fit.temp <- grpreg(cbind(X, Z), Y, lambda = lambda1, 
                           group = c(1:p, rep(0, r)), max.iter = 1e+05, 
                           group.multiplier = w)   
    }
    lambda2 <- list()
    B <- bSpline.ori(U %*% beta0, knots = knots, Boundary.knots = Boundary.knots)$B
    data.temp <- cbind(X[, rep(1:p, each = d)] * B[, rep(1:d, p)], B)
    if (family == "cox") {
        for (lam1 in 1:nlambda) {
            lambda2[[lam1]] <- setup.grlambda(X = data.temp, Y = Y, family = family,
                                              group = c(rep(1:p, each = d), rep(0, d)),
                                              nlambda = nlambda, offset = cbind(X, Z) %*% 
                                                  fit.temp$beta[, lam1], 
                                              lambda.min.ratio = lambda2.min.ratio, 
                                              penalty.factor = c(w, 0))   
        }
    }
    
    if (family == "gaussian") {
        for (lam1 in 1:nlambda) {
            lambda2[[lam1]] <- setup.grlambda(X = data.temp, Y = Y, family = family,
                                              group = c(rep(1:p, each = d), rep(0, d)), 
                                              nlambda = nlambda, offset = cbind(X, Z) %*% 
                                                  fit.temp$beta[-1, lam1], 
                                              lambda.min.ratio = lambda2.min.ratio, 
                                              penalty.factor = c(w, 0))   
        }
    }
    
    # Estimate the unknown parameters
    beta.est <- list()
    coef.est <- list()
    sivcm.likelihood <- matrix(nrow = nlambda, ncol = nlambda)
    sivcm.df <- matrix(nrow = nlambda, ncol = nlambda)
    sivcm.converge <- matrix(nrow = nlambda, ncol = nlambda)
    sivcm.time <- matrix(nrow = nlambda, ncol = nlambda)
    
    beta0.temp <- beta0
    coef0.temp <- coef0
    for (lam1 in nlambda:1) {
        beta.est.temp <- matrix(nrow = nlambda, ncol = q)
        coef.est.temp <- matrix(nrow = nlambda, ncol = (d + 1) * p + d + r + 
                                    as.numeric(family == "gaussian"))
        for (lam2 in nlambda:1) {
            if (lam2 == nlambda) {
                beta0 <- beta0.temp
                coef0 <- coef0.temp
            }
            start.time <- Sys.time()
            est <- try(psivcm(Y, X, U = U, Z = Z, 
                              lambda = c(lambda1[lam1], lambda2[[lam1]][lam2]),
                              family = family, maxit = maxit, tol = 1e-05, 
                              beta0 = beta0, coef0 = coef0, degree = degree,
                              knots = knots, Boundary.knots = Boundary.knots,
                              gm = 1 / w), silent = TRUE)
            end.time <- Sys.time()
            if (class(est) == "try-error")
                next
            if (lam2 == nlambda) {
                beta0.temp <- est$beta
                coef0.temp <- est$coef
            }
            beta0 <- est$beta
            coef0 <- est$coef
            beta.est.temp[lam2, ] <- beta0
            coef.est.temp[lam2, ] <- coef0
            sivcm.likelihood[lam1, lam2] <- est$l
            sivcm.df[lam1, lam2] <- est$df + q
            sivcm.converge[lam1, lam2] <- est$converge
            sivcm.time[lam1, lam2] <- difftime(end.time, start.time, units = "secs")
        }
        beta.est[[lam1]] <- beta.est.temp
        coef.est[[lam1]] <- coef.est.temp
    }
    
    return(list(beta = beta.est, coef = coef.est, likelihood = sivcm.likelihood,
                df = sivcm.df, converge = sivcm.converge, lambda1 = lambda1, 
                lambda2 = lambda2, time = sivcm.time))
}


#-------------------------------------------------------------------------------
#                            Function: tune.control
#-------------------------------------------------------------------------------
# Purpose: A function for model selection by AIC or BIC.
# Arguments:
# family - Response type of "gaussian" or "cox".
# likelihood - Residual sum of squares for "gaussian" or negative log-partial 
#              likelihood for "cox".
# n - Effect sample size.
# df - Effective degrees of freedom.
# tune - Model selection criterion of "bic" (default) or "aic".

tune.control <- function(family, likelihood, n, df, tune = "bic") {
    if (tune == "bic") {
        pen <- log(n) * df 
    }   else if (tune == "aic"){
        pen <- 2 * df
    }
    if (family == "gaussian") {
        ic <- log(likelihood/n) * n + pen
    } else if (family == "cox") {
        ic <- 2 * likelihood + pen
    }
    return(list(lambda.min = which(ic == min(ic, na.rm = TRUE), arr.ind = TRUE),
                ic = ic))
}
