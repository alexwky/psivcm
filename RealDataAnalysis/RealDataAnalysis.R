library(grpreg)
library(psivcm)
library(survival)

source("function.R")

tune <- "bic"
weight <- FALSE


#-------------------------------------------------------------------------------
#                                NSCLC analysis
#-------------------------------------------------------------------------------

# Generate multiple initial values of beta for the NSCLC analysis
if(FALSE) {
    initbeta <- matrix(nrow = 50, ncol = 5)
    for (seed in 1:50) {
        set.seed(seed)
        beta0 <- rnorm(5)
        beta0 <- beta0 / sqrt(sum(beta0 ^ 2))
        initbeta[seed, ] <- beta0
    }
    write.table(initbeta, 
                file = "./RealDataAnalysis/RealDataAnalysisData/RealData-NSCLC-beta0.csv", 
                row.names = FALSE, col.names = FALSE, sep = ",")   
}

# Import data
clinical <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisData/UCSC_NSCLC_FEV1_clinical.csv"))
patient.id <- clinical[, 1]
clinical <- clinical[, -1]
clinical <- apply(clinical, 2, function(x) as.numeric(x))
rownames(clinical) <- patient.id

Y <- clinical[, "FEV1"]
U <- clinical[, c("Age", "PYS", "LUAD", "StageT", "Gender")]
X <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisData/UCSC_NSCLC_FEV1_gene_screen.csv"))
rownames(X) <- X[, 1]
X <- X[, -1]
X <- apply(X, 2, function(x) as.numeric(x))
rownames(X) <- patient.id

n <- length(Y)
p <- ncol(X)
q <- ncol(U)

# Standardize data
X <- apply(X, 2, function(x) (x - mean(x)) / sqrt(sum((x - mean(x)) ^ 2) / n))
U <- apply(U, 2, function(x) (x - mean(x)) / sqrt(sum((x - mean(x)) ^ 2) / n))
Z <- U[, -which(colnames(U) == "StageT")]
r <- ncol(Z)

nknot <- 3
degree <- 2
d <- nknot - 2 + degree
nlambda <- 20

if (weight == TRUE) {
    lambda1.min.ratio <- 5e-07
    lambda2.min.ratio <- 5e-07
} else {
    lambda1.min.ratio <- 0.15
    lambda2.min.ratio <- 0.3
}

# Proposed method
for (initbeta in 1:50) {
    
    if (weight == TRUE) {
        est <- as.matrix(read.csv(
            "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-NSCLC-unweighted.csv",
            sep = ","))
        est <- est[, -1]
        beta0 <- est[initbeta, 1:q]
        coef0 <- est[initbeta, -(1:q)]
        coef0 <- coef0[-length(coef0)]
        gamma.loc <- seq(1, (d + 1) * p, by = d + 1) + 1
        alpha.loc <- (1:((d + 1) * p))[-gamma.loc] + 1
        omega <- 1 / sqrt(rowSums(
            cbind(coef0[gamma.loc], 
                  matrix(coef0[alpha.loc], nrow = p, ncol = d, byrow = T)) ^ 2))
    } else {
        omega <- rep(1, p)
        beta0 <- as.matrix(read.csv(
            "./RealDataAnalysis/RealDataAnalysisData/RealData-NSCLC-beta0.csv",
            header = FALSE))
        beta0 <- beta0[initbeta, ]
        coef0 <- rep(0, (d + 1) * p + d + r + 1)
    }
    
    res <- analysis(Y, X, U, Z, nknot, degree, "gaussian", w = omega, 
                    beta0 = beta0, coef0 = coef0, 
                    lambda1.min.ratio = lambda1.min.ratio, 
                    lambda2.min.ratio = lambda2.min.ratio)
    lambda.loc <- tune.control("gaussian", res$likelihood, n, res$df, tune)$lambda.min
    
    sivcm.beta <- res$beta[[lambda.loc[1]]][lambda.loc[2], ]
    sivcm.coef <- res$coef[[lambda.loc[1]]][lambda.loc[2], ]
    ic <- tune.control("gaussian", res$likelihood, n, res$df, tune)$ic[lambda.loc]
    
    write.table(t(c(initbeta, sivcm.beta, sivcm.coef, ic)), 
                paste("./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-NSCLC-", 
                      ifelse(weight == TRUE, "weighted.csv", "unweighted.csv"),
                      sep = ""),
                sep = ",", row.names = FALSE, 
                col.names = !file.exists(paste(
                    "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-NSCLC-", 
                    ifelse(weight == TRUE, "weighted.csv", "unweighted.csv"),
                    sep = "")), append = T)
}

res <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-NSCLC-weighted.csv",
    row.names = 1))
sivcm.parameter <- res[which.min(res[, ncol(res)]), ]
sivcm.beta <- sivcm.parameter[1:q]
sivcm.coef <- sivcm.parameter[-(1:q)]
sivcm.coef <- sivcm.coef[-length(sivcm.coef)]

# Lasso without interaction terms
if(weight == TRUE){
    coef0 <- as.matrix(read.csv(
        "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-NSCLC-main-unweighted.csv",
        sep = ","))
    omega <- 1 / abs(coef0[(1:p) + 1])
} else {
    omega <- rep(1, p)
}

main.fit <- grpreg(cbind(X, U), Y, group = c(1:p, rep(0, q)), family = "gaussian",
                   lambda = setup.grlambda(cbind(X, U), Y, family = "gaussian",
                                           group = c(1:p, rep(0, q)), nlambda = nlambda,
                                           lambda.min.ratio = lambda1.min.ratio,
                                           penalty.factor = c(omega, 0)),
                   group.multiplier = omega)
main.likelihood <- apply(main.fit$beta, 2, function(x)
    likelihood(Y, colSums(t(cbind(1, X, U)) * x), family = "gaussian"))
main.tune <- tune.control(family = "gaussian", main.likelihood, n, main.fit$df, tune)
main.coef <- main.fit$beta[,main.tune$lambda.min]
write.table(t(main.coef), 
            paste("./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-NSCLC-main-", 
                  ifelse(weight == TRUE, "weighted.csv", "unweighted.csv"),
                  sep = ""), sep = ",", row.names = FALSE)


#-------------------------------------------------------------------------------
#                                 LGG analysis
#-------------------------------------------------------------------------------

# Generate multiple initial values of beta for the LGG analysis
if(FALSE) {
    initbeta <- matrix(nrow = 50, ncol = 7)
    for (seed in 1:50) {
        set.seed(seed)
        beta0 <- rnorm(7)
        beta0 <- beta0 / sqrt(sum(beta0 ^ 2))
        initbeta[seed, ] <- beta0
    }
    write.table(initbeta, 
                file = "./RealDataAnalysis/RealDataAnalysisData/RealData-LGG-beta0.csv", 
                row.names = FALSE, col.names = FALSE, sep = ",")
}

# Import data
Y <- as.matrix(read.csv("./RealDataAnalysis/RealDataAnalysisData/UCSC_LGG_OS_survival.csv"))
patient.id <- Y[, 1]
Y <- Surv(as.numeric(Y[, 2]), as.numeric(Y[, 3]))
names(Y) <- patient.id

X <- as.matrix(read.csv("./RealDataAnalysis/RealDataAnalysisData/UCSC_LGG_OS_protein.csv"))
rownames(X) <- X[, 1]
X <- X[, -1]
X <- apply(X, 2, function(x) as.numeric(x))
rownames(X) <- patient.id

U <- as.matrix(read.csv("./RealDataAnalysis/RealDataAnalysisData/UCSC_LGG_OS_gene.csv"))
rownames(U) <- U[, 1]
U <- U[, -1]
U <- apply(U, 2, function(x) as.numeric(x))
rownames(U) <- patient.id

Z <- as.matrix(read.csv("./RealDataAnalysis/RealDataAnalysisData/UCSC_LGG_OS_clinical.csv"))
rownames(Z) <- Z[, 1]
Z <- Z[, -1]
Z <- apply(Z, 2, function(x) as.numeric(x))
rownames(Z) <- patient.id

npc <- 7
pc <- prcomp(U, rank. = npc)
U <- pc$x

n <- length(Y)
p <- ncol(X)
q <- ncol(U)

# Standardize data
X <- apply(X, 2, function(x) (x - mean(x)) / sqrt(sum((x - mean(x)) ^ 2) / n))
U <- apply(U, 2, function(x) (x - mean(x)) / sqrt(sum((x - mean(x)) ^ 2) / n))
Z <- apply(Z, 2, function(x) (x - mean(x)) / sqrt(sum((x - mean(x)) ^ 2) / n))
Z <- cbind(Z, U[, -q])
r <- ncol(Z)

ord <- order(Y[, 1], -Y[, 2])
Y <- Y[ord, ]
X <- X[ord, ]
U <- U[ord, ]
Z <- Z[ord, ]

nknot <- 3
degree <- 2
d <- nknot - 2 + degree
nlambda <- 20

if (weight == TRUE) {
    lambda1.min.ratio <- 5e-07
    lambda2.min.ratio <- 5e-07
} else {
    lambda1.min.ratio <- 0.15
    lambda2.min.ratio <- 0.15
}

# Proposed method
for (initbeta in 1:50) {
    
    if (weight == TRUE) {
        est <- as.matrix(read.csv(
            "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-LGG-unweighted.csv",
            sep = ","))
        est <- est[, -1]
        beta0 <- est[initbeta, 1:q]
        coef0 <- est[initbeta, -(1:q)]
        coef0 <- coef0[- length(coef0)]
        gamma.loc <- seq(1, (d + 1) * p, by = d + 1)
        alpha.loc <- (1:((d + 1) * p))[-gamma.loc]
        omega <- 1 / sqrt(rowSums(
            cbind(coef0[gamma.loc], 
                  matrix(coef0[alpha.loc], nrow = p, ncol = d, byrow = T)) ^ 2))
        omega[omega == Inf] <- 1 / 1e-300
    } else {
        omega <- rep(1, p)
        beta0 <- as.matrix(read.csv(
            "./RealDataAnalysis/RealDataAnalysisData/RealData-LGG-beta0.csv",
            header = FALSE))
        beta0 <- beta0[initbeta, ]
        coef0 <- rep(0, (d + 1) * p + d + r)
    }
    
    res <- analysis(Y, X, U, Z, nknot, degree, "cox", w = omega, beta0 = beta0,
                    coef0 = coef0, lambda1.min.ratio = lambda1.min.ratio, 
                    lambda2.min.ratio = lambda2.min.ratio, maxit = 500)
    lambda.loc <- tune.control("cox", res$likelihood, sum(Y[, 2]), res$df, tune)$lambda.min
    
    sivcm.beta <- res$beta[[lambda.loc[1]]][lambda.loc[2], ]
    sivcm.coef <- res$coef[[lambda.loc[1]]][lambda.loc[2], ]
    ic <- tune.control("cox", res$likelihood, sum(Y[, 2]), res$df, tune)$ic[lambda.loc]
    
    write.table(t(c(initbeta, sivcm.beta, sivcm.coef, ic)),
                paste("./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-LGG-",
                      ifelse(weight == TRUE, "weighted.csv", "unweighted.csv"),
                      sep = ""),
                sep = ",", row.names = FALSE, 
                col.names = !file.exists(paste(
                    "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-LGG-",
                    ifelse(weight == TRUE, "weighted.csv", "unweighted.csv"),
                    sep = "")), append = T)
}

res <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-LGG-weighted.csv",
    row.names = 1))
sivcm.parameter <- res[which.min(res[, ncol(res)]), ]
sivcm.beta <- sivcm.parameter[1:q]
sivcm.coef <- sivcm.parameter[-(1:q)]
sivcm.coef <- sivcm.coef[-length(sivcm.coef)]

# Lasso without interaction terms
if(weight == TRUE){
    coef0 <- as.matrix(read.csv(
        "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-LGG-main-unweighted.csv",
        sep = ","))
    omega <- 1 / abs(coef0[(1:p)])
    omega[omega == Inf] <- 1 / 1e-300
} else {
    omega <- rep(1, p)
}

main.fit <- grpsurv(cbind(X, U), Y, group = c(1:p, rep(0, q)), family = "cox",
                    lambda = setup.grlambda(cbind(X, U), Y, family = "cox",
                                            group = c(1:p, rep(0, q)), nlambda = nlambda,
                                            lambda.min.ratio = lambda1.min.ratio,
                                            penalty.factor = c(omega, 0)),
                    group.multiplier = omega)
main.likelihood <- apply(main.fit$beta, 2, function(x)
    likelihood(Y, colSums(t(cbind(X, U)) * x), family = "cox"))
main.tune <- tune.control(family = "cox", main.likelihood, sum(Y[,2]), main.fit$df,
                          tune)
main.coef <- main.fit$beta[,main.tune$lambda.min]
write.table(t(main.coef), 
            paste("./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-LGG-main-", 
                  ifelse(weight == TRUE, "weighted.csv", "unweighted.csv"),
                  sep = ""), sep = ",", row.names = FALSE)
