#-------------------------------------------------------------------------------
#                                    Setting
#-------------------------------------------------------------------------------

method <- vector(length = 3)
names(method) <- c("SIVC", "MAIN", "INT")
method["SIVC"] <- TRUE
method["MAIN"] <- TRUE # FALSE for model 3
method["INT"] <- TRUE # FALSE for model 3

# --------------------------------- User inputs --------------------------------
p <- 20 # number of covariates in X: 20, 50, 100, and 300
family <- "gaussian" # response type: "gaussian" and "cox"
model <- 1 # simulation setting: 1, 2, 3, and 4
nonlinear.main <- FALSE # logical argument for fitting nonlinear main effect
tune <- "bic" # model selection criterion: "bic", "aic"
initbeta <- 1 # set of initial value of beta: 1, 2, 3, 4, and 5
weight <- FALSE # logical argument for the weighted method
# (users should first run the analysis for weight = FALSE)
# ------------------------------------------------------------------------------

n <- 500
q <- 4
beta <- c(0.4, -0.4, 0.2, 0.8)
psi <- c(0.2, -0.2, 0.5, -0.5)
{
    g1 <- function(z) -0.5 + 0 * z
    g2 <- function(z) -0.4 + 0 * z
    g3 <- function(z) -0.3 + 0 * z
    g4 <- function(z) -0.2 + 0 * z
    g5 <- function(z) -0.1 + 0 * z
    g6 <- function(z) 0.1 + 0 * z
    g7 <- function(z) 0.2 + 0 * z
    g8 <- function(z) 0.3 + 0 * z
    g9 <- function(z) 0.4 + 0 * z
    g10 <- function(z) 0.5 + 0 * z
    if (model == 2) {
        g11 <- function(z) -0.5 + 0 * z
        g12 <- function(z) -0.4 + 0 * z
        g13 <- function(z) -0.3 + 0 * z
        g14 <- function(z) -0.2 + 0 * z
        g15 <- function(z) -0.1 + 0 * z
        g16 <- function(z) 0.1 + 0 * z
        g17 <- function(z) 0.2 + 0 * z
        g18 <- function(z) 0.3 + 0 * z
        g19 <- function(z) 0.4 + 0 * z
        g20 <- function(z) 0.5 + 0 * z
    } else {
        g11 <- function(z) 0.25 * z + 0 * z ^ 2
        g12 <- function(z) -0.2 - 0.3 * z + 0 * z ^ 2
        g13 <- function(z) -0.25 * z ^ 2
        g14 <- function(z) -0.2 + 0.25 * z ^ 2
        g15 <- function(z) 0.4 * z - 0.1 * z ^ 2
        g16 <- function(z) -0.3 + 0.3 * z + 0.2 * z ^ 2
        g17 <- function(z) 0.1 * z ^ 3
        g18 <- function(z) -0.2 - 0.15 * z ^ 3
        g19 <- function(z) 1 / (1 + exp(- z * 5))
        g20 <- function(z) 1 / (1 + exp(z * 4)) - 0.3
    }
}
nknot <- 3
degree <- 2
d <- nknot - 2 + degree
nlambda <- 20
if (weight) {
    lambda1.min.ratio <- 0.001
    lambda2.min.ratio <- 0.001
} else {
    if (p == 20) {
        lambda1.min.ratio <- 0.001
        lambda2.min.ratio <- 0.001
    } else {
        lambda1.min.ratio <- 0.05
        lambda2.min.ratio <- 0.05
    }
}


#-------------------------------------------------------------------------------
#                                 Run analyses                                 
#-------------------------------------------------------------------------------

library(psivcm)
library(grpreg)
library(survival)
library(survC1)

source("function.R")

for (seed in 1:100) {
    
    data.train <- as.matrix(read.csv(
        paste("./Simulation/SimulationData/SimulationData-p", p, "-", seed,
              ".csv", sep = "")))
    X <- data.train[, which(sapply(1:ncol(data.train), function(x) 
        unlist(strsplit(colnames(data.train)[x], ""))[1]) == "X")]
    U <- data.train[, which(sapply(1:ncol(data.train), function(x) 
        unlist(strsplit(colnames(data.train)[x], ""))[1]) == "U")]
    if (model == 3) {
        if (nonlinear.main == TRUE) {
            Z <- cbind(U ^ 2, U)
            Z <- Z[, -ncol(Z)]
        } else Z <- U[, -q]
    } else Z <- U[, -q]
    r <- ncol(Z)
    if (model == 1) {
        if (family == "cox") Y <- Surv(data.train[, "time"], data.train[, "status"])
        if (family == "gaussian")  Y <- data.train[, "Y"]
    } else {
        response.train <- as.matrix(read.csv(
            paste("./Simulation/SimulationData/AdditionalSimulationData-setting",
                  model, "-p", p, "-", seed, ".csv", sep = "")))
        if (family == "cox") Y <- Surv(response.train[, "time"], response.train[, "status"])
        if (family == "gaussian") Y <- response.train[, "Y"]
        if (model == 2) {
            Z <- matrix(nrow = n, ncol = 0)
            r <- ncol(Z)
        }
    }
    
    if (family == "cox") {
        ord <- order(Y[, 1], -Y[, 2])
        Y <- Y[ord, ]
        X <- X[ord, ]
        U <- U[ord, ]
        Z <- Z[ord, ]
    }
    
    data.val <- as.matrix(read.csv(
        paste("./Simulation/SimulationData/SimulationData-p", p, "-", 0, ".csv",
              sep = "")))
    X.val <- data.val[, which(sapply(1:ncol(data.val), function(x) 
        unlist(strsplit(colnames(data.val)[x], ""))[1]) == "X")]
    U.val <- data.val[, which(sapply(1:ncol(data.val), function(x) 
        unlist(strsplit(colnames(data.val)[x], ""))[1]) == "U")]
    if (model == 3) {
        if (nonlinear.main == TRUE) {
            Z.val <- cbind(U.val ^ 2, U.val)
            Z.val <- Z.val[, -ncol(Z.val)]
        } else Z.val <- U.val[, -q]
    } else Z.val <- U.val[, -q]
    
    g.val <- matrix(0, nrow = 5000, ncol = p)
    g.val[, 1:20] <- sapply(1:20, function(i) 
        do.call(paste("g", i, sep = ""), list(U.val %*% beta)))
    g.val <- cbind(g.val, matrix(psi, nrow = 5000, ncol = q, byrow = TRUE))
    if (model == 3) {
        U.val.trans <- U.val
        U.val.trans[, 1] <- 0.5 * U.val[, 1] ^ 2 + 0.8 * U.val[, 1]
        U.val.trans[, 2] <- 0.5 * U.val[, 2] ^ 2 + 0.8 * U.val[, 2]
        U.val.trans[, 3] <- -0.3 * U.val[, 3] ^ 2 + 0.8 * U.val[, 3]
        U.val.trans[, 4] <- -0.3 * U.val[, 4] ^ 2 + 0.8 * U.val[, 4]
        eta <- rowSums(g.val * cbind(X.val, U.val.trans))
    } else eta <- rowSums(g.val * cbind(X.val, U.val))
    
    if (model == 1) {
        if (family == "cox") Y.val <- Surv(data.val[, "time"], data.val[, "status"])
        if (family == "gaussian") Y.val <- data.val[, "Y"]
    } else {
        response.val <- as.matrix(read.csv(
            paste("./Simulation/SimulationData/AdditionalSimulationData-setting",
                  model, "-p", p, "-", 0, ".csv", sep = "")))
        if (family == "cox") Y.val <- Surv(response.val[, "time"], response.val[, "status"])
        if (family == "gaussian") Y.val <- response.val[, "Y"]
        if (model == 2) Z.val <- matrix(nrow = 5000, ncol = 0)
    }
    
    # Proposed method
    if (method["SIVC"]) {
        if (weight == TRUE) {
            est <- as.matrix(read.csv(
                paste("./Simulation/SimulationResults/Estimates-setting", model, 
                      ifelse(nonlinear.main == TRUE, "-", ifelse(model == 3, "-mis-", "-")),
                      family, "-p", p, "-beta", initbeta,
                      ifelse(tune == "aic", "-unweighted-aic.csv", "-unweighted.csv"), 
                      sep = "")))
            beta0 <- est[seed, 1:q]
            coef0 <- est[seed, -(1:q)]
            gamma.loc <- seq(1, (d + 1) * p, by = d + 1)
            alpha.loc <- (1:((d + 1) * p))[-gamma.loc]
            if (family == "gaussian") {
                gamma.loc <- gamma.loc + 1
                alpha.loc <- alpha.loc + 1
            }
            omega <- 1 / sqrt(rowSums(cbind(coef0[gamma.loc], 
                                            matrix(coef0[alpha.loc], nrow = p, 
                                                   ncol = d, byrow = T)) ^ 2))
            if (family == "cox") omega[omega == Inf] <- 1 / 1e-300
        } else {
            omega <- rep(1, p)
            beta0 <- as.matrix(read.csv(
                "./Simulation/SimulationData/Simulation-beta0.csv", header = FALSE))
            beta0 <- beta0[initbeta, ]
            coef0 <- rep(0, (d + 1) * p + d + r + as.numeric(family == "gaussian"))
        }
        
        res <- analysis(Y, X, U, Z, nknot, degree, family, w = omega, beta0 = beta0,
                        coef0 = coef0, lambda1.min.ratio = lambda1.min.ratio, 
                        lambda2.min.ratio = lambda2.min.ratio)
        lambda.loc <- tune.control(family, res$likelihood, 
                                   ifelse(family == "gaussian", n, sum(Y[, 2])), 
                                   res$df, tune)$lambda.min
        sivcm.beta <- res$beta[[lambda.loc[1]]][lambda.loc[2], ]
        sivcm.coef <- res$coef[[lambda.loc[1]]][lambda.loc[2], ]
        
        # Evaluate the model performance
        sivcm.result <- vector(length = ifelse(model == 2, 12, 15))
        if (model == 2) {
            names(sivcm.result) <- c("seed", "sen", "fdr", "card", "card.main",
                                     "card.int", "mse", "c-index", "lambda1",
                                     "lambda2", "converge", "ic")
        } else {
            names(sivcm.result) <- c("seed", "sen", "sen.int", "fdr", "fdr.int",
                                     "card", "card.main", "card.int", "aip",
                                     "mse", "c-index", "lambda1", "lambda2",
                                     "converge", "ic")
        }
        sivcm.result["aip"] <- abs(beta %*% sivcm.beta)
        sel.true <- 1:20
        sel.int.true <- 11:20
            
        nsel.true <- (1:p)[-sel.true]
        nsel.int.true <- (1:p)[-sel.int.true]
        sivcm.result["seed"] <- seed
        sivcm.coef.mat <- matrix(sivcm.coef[(1:(p * (d + 1))) + 
                                                as.numeric(family == "gaussian")], 
                                 ncol = d + 1, byrow = TRUE)
        sel.sivcm <- which(apply(sivcm.coef.mat, 1, function(x) !all(x == 0)))
        sel.int.sivcm <- which(apply(sivcm.coef.mat[, -1], 1, function(x) any(x != 0)))
        sivcm.result["sen"] <- length(intersect(sel.sivcm, sel.true)) / length(sel.true)
        if (length(intersect(nsel.true, sel.sivcm)) == 0) {
            sivcm.result["fdr"] <- 0
        } else {
            sivcm.result["fdr"] <- length(intersect(nsel.true, sel.sivcm)) / 
                length(sel.sivcm)
        }
        sivcm.result["sen.int"] <- length(intersect(sel.int.sivcm, sel.int.true)) / 
            length(sel.int.true)
        if (length(intersect(nsel.int.true, sel.int.sivcm)) == 0) {
            sivcm.result["fdr.int"] <- 0
        } else {
            sivcm.result["fdr.int"] <- length(intersect(nsel.int.true, sel.int.sivcm)) / 
                length(sel.int.sivcm)
        }
        sivcm.result[c("lambda1", "lambda2")] <- c(res$lambda1[lambda.loc[1]], 
                                                   res$lambda2[[lambda.loc[1]]][lambda.loc[2]])
        sivcm.result["card"] <- length(sel.sivcm)
        sivcm.result["card.main"] <- sum(sivcm.coef.mat[, 1] != 0)
        sivcm.result["card.int"] <- length(sel.int.sivcm)
        sivcm.result["converge"] <- res$converge[lambda.loc[1], lambda.loc[2]]
        sivcm.result["ic"] <- tune.control(family, res$likelihood, 
                                           ifelse(family == "gaussian", n, sum(Y[, 2])),
                                           res$df, tune)$ic[lambda.loc]
        knots <- qnorm(seq(0, 1, l = nknot)[-c(1, nknot)])
        Boundary.knots <- c(-max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))), 
                            max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))))
        B <- bSpline.ori(U.val %*% sivcm.beta, knots = knots, 
                         Boundary.knots = Boundary.knots)$B
        sivcm.g.val <- get.g(X.val, U.val, Z.val, B, sivcm.beta, sivcm.coef, 
                             family = family)
        if (family == "cox") {
            sivcm.result["mse"] <- mean((eta - rowSums(
                sivcm.g.val * cbind(X.val, Z.val, B))) ^ 2)
        }
        if (family == "gaussian") {
            sivcm.result["mse"] <- mean((eta - rowSums(
                sivcm.g.val * cbind(1, X.val, Z.val, B))) ^ 2)   
        }
        if (family == "cox") {
            sivcm.result["c-index"] <- 
                Est.Cval(mydata = cbind(Y.val[, 1] * 100, Y.val[, 2], 
                                        sapply(1:5000, function(i) sivcm.g.val[i, ] %*%
                                                   c(X.val[i, ], Z.val[i, ], B[i, ]))), 
                         tau = max(Y.val[, 1]) * 100, nofit = TRUE)$Dhat
        } else if (family == "gaussian") sivcm.result["c-index"] <- NA
        sivcm.result <- sivcm.result[1:ifelse(model == 2, 12, length(sivcm.result))]
        
        write.table(t(sivcm.result), paste(
            "./Simulation/SimulationResults/SimulationResults-setting", model,
            ifelse(nonlinear.main == TRUE, "-", ifelse(model == 3, "-mis-", "-")), 
            family, "-p", p, "-beta", initbeta, "-", 
            ifelse(weight == TRUE, "weighted", "unweighted"), 
            ifelse(tune == "aic", "-aic.csv", ".csv"), sep = ""),
            sep = ",", row.names = FALSE,
            col.names = !file.exists(paste(
                "./Simulation/SimulationResults/SimulationResults-setting", model, 
                ifelse(nonlinear.main == TRUE, "-", ifelse(model == 3, "-mis-", "-")),
                family, "-p", p, "-beta", initbeta, "-", 
                ifelse(weight == TRUE, "weighted", "unweighted"),
                ifelse(tune == "aic", "-aic.csv", ".csv"), sep = "")), append = T)
        write.table(t(c(sivcm.beta, sivcm.coef)), paste(
            "./Simulation/SimulationResults/Estimates-setting", model, 
            ifelse(nonlinear.main == TRUE, "-", ifelse(model == 3, "-mis-", "-")),
            family, "-p", p, "-beta", initbeta, "-", 
            ifelse(weight == TRUE, "weighted", "unweighted"),
            ifelse(tune == "aic", "-aic.csv", ".csv"), sep = ""), 
            sep = ",", row.names = FALSE,
            col.names = !file.exists(paste(
                "./Simulation/SimulationResults/Estimates-setting", model,
                ifelse(nonlinear.main == TRUE, "-", ifelse(model == 3, "-mis-", "-")),
                family, "-p", p, "-beta", initbeta, "-",
                ifelse(weight == TRUE, "weighted", "unweighted"),
                ifelse(tune == "aic", "-aic.csv", ".csv"), sep = "")), append = T)
    }
    
    # Lasso without interaction terms
    if (method["MAIN"]) {
        if (weight == TRUE) {
            est <- as.matrix(read.csv(paste(
                "./Simulation/SimulationResults/Estimates-setting",
                model, "-", family, "-p", p, 
                ifelse(tune == "aic", "-main-unweighted-aic.csv", "-main-unweighted.csv"),
                sep = "")), header = FALSE)
            if (family == "cox") {
                omega <- 1 / abs(est[seed, 1:p])
                omega[omega == Inf] <- 1 / 1e-300
            }
            if (family == "gaussian") omega <- 1 / abs(est[seed, (1:p) + 1])
        } else {
            omega <- rep(1, p)
        }
        
        if (family == "cox") {
            main.fit <- grpsurv(cbind(X, U), Y, group = c(1:p, rep(0, q)), 
                                lambda = setup.grlambda(cbind(X, U), Y, family = family,
                                                        group = c(1:p, rep(0, q)),
                                                        nlambda = nlambda,
                                                        lambda.min.ratio = lambda1.min.ratio,
                                                        penalty.factor = c(omega, 0)),
                                group.multiplier = omega)
            main.likelihood <- apply(main.fit$beta, 2, function(x) 
                likelihood(Y, colSums(t(cbind(X, U)) * x), family))
        }
        if (family == "gaussian") {
            main.fit <- grpreg(cbind(X, U), Y, group = c(1:p, rep(0, q)), family = family,
                               lambda = setup.grlambda(cbind(X, U), Y, family = family,
                                                       group = c(1:p, rep(0, q)),
                                                       nlambda = nlambda,
                                                       lambda.min.ratio = lambda1.min.ratio,
                                                       penalty.factor = c(omega, 0)),
                               group.multiplier = omega)
            main.likelihood <- apply(main.fit$beta, 2, function(x) 
                likelihood(Y, colSums(t(cbind(1, X, U)) * x), family))  
        }
        main.tune <- tune.control(family, main.likelihood, 
                                  ifelse(family == "gaussian", n, sum(Y[, 2])), 
                                  main.fit$df, tune)
        main.coef <- main.fit$beta[, main.tune$lambda.min]
        main.result <- vector(length = 7)
        names(main.result) <- c("seed", "sen", "fdr", "card", "mse", "c-index", "lambda")
        if (model == 4) sel.true <- 1:2 else sel.true <- 1:20
        nsel.true <- (1:p)[-sel.true]
        main.result["seed"] <- seed
        if (family == "cox") {
            sel.main <- which(main.fit$beta[1:p, main.tune$lambda.min] != 0)  
        }
        if (family == "gaussian") {
            sel.main <- which(main.fit$beta[(1:p) + 1, main.tune$lambda.min] != 0)   
        }
        main.result["sen"] <- length(intersect(sel.main, sel.true)) / length(sel.true)
        if (length(intersect(nsel.true, sel.main)) == 0) {
            main.result["fdr"] <- 0
        } else {
            main.result["fdr"] <- length(intersect(nsel.true, sel.main)) / length(sel.main)
        }
        main.result["lambda"] <- main.fit$lambda[main.tune$lambda.min]
        main.result["card"] <- length(sel.main)
        if (family == "cox") {
            main.g <- matrix(main.fit$beta[, main.tune$lambda.min], nrow = 5000,
                             ncol = p + q, byrow = TRUE)
            main.result["mse"] <- mean((eta - rowSums(main.g * cbind(X.val, U.val))) ^ 2)
            main.result["c-index"] <- Est.Cval(
                mydata = cbind(Y.val[, 1] * 100, Y.val[, 2], sapply(1:5000, function(i) 
                    main.g[i, ] %*% c(X.val[i, ], U.val[i, ]))),
                tau = max(Y.val[, 1]) * 100, nofit = TRUE)$Dhat
        }
        if (family == "gaussian") {
            main.g <- matrix(main.fit$beta[, main.tune$lambda.min], nrow = 5000,
                             ncol = p + q + 1, byrow = TRUE)
            main.result["mse"] <- mean((eta - rowSums(main.g * cbind(1, X.val, U.val))) ^ 2)
            main.result["c-index"] <- NA
        }
        
        write.table(t(main.result), paste(
            "./Simulation/SimulationResults/SimulationResults-setting",
            model, "-", family, "-p", p, "-main-",
            ifelse(weight == TRUE, "weighted", "unweighted"),
            ifelse(tune == "aic", "-aic.csv", ".csv"), sep = ""),
            sep = ",", row.names = FALSE,
            col.names = !file.exists(paste(
                "./Simulation/SimulationResults/SimulationResults-setting",
                model, "-", family, "-p", p, "-main-",
                ifelse(weight == TRUE, "weighted", "unweighted"),
                ifelse(tune == "aic", "-aic.csv", ".csv"), sep = "")), append = T)
        write.table(t(main.coef), paste(
            "./Simulation/SimulationResults/Estimates-setting",
            model, "-", family, "-p", p, "-main-",
            ifelse(weight == TRUE, "weighted", "unweighted"),
            ifelse(tune == "aic", "-aic.csv", ".csv"), sep = ""),
            sep = ",", row.names = FALSE,
            col.names = !file.exists(paste(
                "./Simulation/SimulationResults/Estimates-setting",
                model, "-", family, "-p", p, "-main-",
                ifelse(weight == TRUE, "weighted", "unweighted"),
                ifelse(tune == "aic", "-aic.csv", ".csv"), sep = "")), append = T)
    }
    
    # Lasso with interaction terms
    if (method["INT"]) {
        if (weight == TRUE) {
            est <- as.matrix(read.csv(paste(
                "./Simulation/SimulationResults/Estimates-setting",
                model, "-", family, "-p", p,
                ifelse(tune == "aic", "-int-unweighted-aic.csv", "-int-unweighted.csv"), 
                sep = "")), header = FALSE)
            if (family == "cox") {
                omega <- 1 / abs(est[seed, (1:(p * (q + 1)))])
                omega[omega == Inf] <- 1 / 1e-300
            }
            if (family == "gaussian") omega <- 1 / abs(est[seed, (1:(p * (q + 1))) + 1])
        } else {
            omega <- rep(1, p * (q + 1))
        }
        
        if (family == "cox") {
            int.fit <- grpsurv(cbind(
                X[, rep(1:p, each = q + 1)] *
                    matrix(cbind(rep(1, n), U), nrow = n, ncol = p * (q + 1)), U),
                Y, group = c(1:(p * (q + 1)), rep(0, q)), 
                lambda = setup.grlambda(cbind(
                    X[, rep(1:p, each = q + 1)] *
                        matrix(cbind(rep(1, n), U), nrow = n, ncol = p * (q + 1)),
                    U), Y, family = family, group = c(1:(p * (q + 1)), rep(0, q)),
                    nlambda = nlambda, lambda.min.ratio = lambda1.min.ratio,
                    penalty.factor = c(omega, 0)), group.multiplier = omega)
            int.likelihood <- apply(int.fit$beta, 2, function(x) 
                likelihood(Y, colSums(t(cbind(X[, rep(1:p, each = q + 1)] * 
                                                  matrix(cbind(rep(1, n), U), nrow = n, 
                                                         ncol = p * (q + 1)), U)) * x), family))
        }
        
        if (family == "gaussian") {
            int.fit <- grpreg(cbind(X[, rep(1:p, each = q + 1)] * 
                                        matrix(cbind(rep(1, n), U), nrow = n, ncol = p * (q + 1)),
                                    U), Y, family = family, group = c(1:(p * (q + 1)), rep(0, q)),
                              lambda = setup.grlambda(cbind(
                                  X[, rep(1:p, each = q + 1)] *
                                      matrix(cbind(rep(1, n), U), nrow = n, ncol = p * (q + 1)),
                                  U), Y, family = family, group = c(1:(p * (q + 1)), rep(0, q)),
                                  nlambda = nlambda, lambda.min.ratio = lambda1.min.ratio,
                                  penalty.factor = c(omega, 0)), group.multiplier = omega)
            int.likelihood <- apply(int.fit$beta, 2, function(x) 
                likelihood(Y, colSums(t(cbind(1, X[, rep(1:p, each = q + 1)] * 
                                                  matrix(cbind(rep(1, n), U), nrow = n, 
                                                         ncol = p * (q + 1)), U)) * x), family))
        }
        
        int.tune <- tune.control(family, int.likelihood, 
                                 ifelse(family == "gaussian", n, sum(Y[, 2])),
                                 int.fit$df, tune)
        int.coef <- int.fit$beta[, int.tune$lambda.min]
        if (family == "cox") {
            int.coef.mat <- matrix(int.coef[1:(p * (q + 1))], ncol = q + 1, byrow = TRUE)
        }
        if (family == "gaussian") {
            int.coef.mat <- matrix(int.coef[(1:(p * (q + 1))) + 1], ncol = q + 1, byrow = TRUE)
        }
        int.result <- vector(length = ifelse(model == 2, 9, 11))
        if (model == 2) {
            names(int.result) <- c("seed", "sen", "fdr", "card", "card.main",
                                   "card.int", "mse", "c-index", "lambda")
        } else {
            names(int.result) <- c("seed", "sen", "sen.int", "fdr", "fdr.int",
                                   "card", "card.main", "card.int", "mse",
                                   "c-index", "lambda")
        }
        sel.true <- 1:20
        sel.int.true <- 11:20
        nsel.true <- (1:p)[-sel.true]
        nsel.int.true <- (1:p)[-sel.int.true]
        sel.int <- which(apply(int.coef.mat, 1, function(x) any(x != 0)))
        sel.int.int <- which(apply(int.coef.mat[, -1], 1, function(x) any(x != 0)))
        int.result["seed"] <- seed
        int.result["sen"] <- length(intersect(sel.int, sel.true)) / length(sel.true)
        if (length(intersect(nsel.true, sel.int)) == 0) {
            int.result["fdr"] <- 0
        } else {
            int.result["fdr"] <- length(intersect(nsel.true, sel.int)) / length(sel.int)
        }
        int.result["sen.int"] <- length(intersect(sel.int.int, sel.int.true)) / 
            length(sel.int.true)
        if (length(intersect(nsel.int.true, sel.int.int)) == 0) {
            int.result["fdr.int"] <- 0
        } else {
            int.result["fdr.int"] <- length(intersect(nsel.int.true, sel.int.int)) / 
                length(sel.int.int)
        }
        int.result["lambda"] <- int.fit$lambda[int.tune$lambda.min]
        int.result["card"] <- length(sel.int)
        int.result["card.main"] <- sum(int.coef.mat[, 1] != 0)
        int.result["card.int"] <- length(sel.int.int)
        if (family == "cox") {
            int.g <- matrix(nrow = 5000, ncol = p + q)
            for (i in 1:5000) {
                temp.mat <- matrix(int.coef[1:(p * (q + 1))], nrow = p, ncol = q +
                                       1, byrow = TRUE)
                temp.mat[, 1] <- temp.mat[, 1]
                temp.mat[, 2] <- temp.mat[, 2] * U.val[i, 1]
                temp.mat[, 3] <- temp.mat[, 3] * U.val[i, 2]
                temp.mat[, 4] <- temp.mat[, 4] * U.val[i, 3]
                temp.mat[, 5] <- temp.mat[, 5] * U.val[i, 4]
                int.g[i, ] <- c(rowSums(temp.mat), tail(int.coef, q))
            }
            int.result["mse"] <- mean((eta - rowSums(int.g * cbind(X.val, U.val))) ^ 2)
            int.result["c-index"] <- Est.Cval(
                mydata = cbind(Y.val[, 1] * 100, Y.val[, 2],
                               sapply(1:5000, function(i) int.g[i, ] %*% 
                                          c(X.val[i, ], U.val[i, ]))), 
                tau = max(Y.val[, 1]) * 100, nofit = TRUE)$Dhat
        } else if (family == "gaussian") {
            int.g <- matrix(nrow = 5000, ncol = p + q + 1)
            for (i in 1:5000) {
                temp.mat <- matrix(int.coef[(1:(p * (q + 1))) + 1], 
                                   nrow = p, ncol = q + 1, byrow = TRUE)
                temp.mat[, 1] <- temp.mat[, 1]
                temp.mat[, 2] <- temp.mat[, 2] * U.val[i, 1]
                temp.mat[, 3] <- temp.mat[, 3] * U.val[i, 2]
                temp.mat[, 4] <- temp.mat[, 4] * U.val[i, 3]
                temp.mat[, 5] <- temp.mat[, 5] * U.val[i, 4]
                int.g[i, ] <- c(int.coef[1], rowSums(temp.mat), tail(int.coef, q))
            }
            int.result["mse"] <- mean((eta - rowSums(int.g * cbind(1, X.val, U.val))) ^ 2)
            int.result["c-index"] <- NA
        }
        int.result <- int.result[1:ifelse(model == 2, 9, length(int.result))]
        
        write.table(t(int.result), paste(
            "./Simulation/SimulationResults/SimulationResults-setting",
            model, "-", family, "-p", p, "-int-",
            ifelse(weight == TRUE, "weighted", "unweighted"),
            ifelse(tune == "aic", "-aic.csv", ".csv"), sep = ""),
            sep = ",", row.names = FALSE,
            col.names = !file.exists(paste(
                "./Simulation/SimulationResults/SimulationResults-setting",
                model, "-", family, "-p", p, "-int-",
                ifelse(weight == TRUE, "weighted", "unweighted"),
                ifelse(tune == "aic", "-aic.csv", ".csv"), sep = "")), append = T)
        write.table(t(int.coef), paste(
            "./Simulation/SimulationResults/Estimates-setting",
            model, "-", family, "-p", p, "-int-",
            ifelse(weight == TRUE, "weighted", "unweighted"),
            ifelse(tune == "aic", "-aic.csv", ".csv"), sep = ""),
            sep = ",", row.names = FALSE,
            col.names = !file.exists(paste(
                "./Simulation/SimulationResults/Estimates-setting",
                model, "-", family, "-p", p, "-int-",
                ifelse(weight == TRUE, "weighted", "unweighted"),
                ifelse(tune == "aic", "-aic.csv", ".csv"), sep = "")), append = T)
    }
}
