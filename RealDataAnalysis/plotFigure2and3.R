library(psivcm)

#-------------------------------------------------------------------------------
#                                    Figure 2
#-------------------------------------------------------------------------------

# Import data
clinical <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisData/UCSC_NSCLC_FEV1_clinical.csv"))
patient.id <- clinical[, 1]
clinical <- clinical[, -1]
clinical <- apply(clinical, 2, function(x) as.numeric(x))
rownames(clinical) <- patient.id

Y <- clinical[, "FEV1"]
U <- clinical[, c("Age", "PYS", "LUAD", "Gender", "StageT")]
X <- as.matrix(read.csv("./RealDataAnalysis/RealDataAnalysisData/UCSC_NSCLC_FEV1_gene_screen.csv"))
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
Z <- U[, -q]
r <- ncol(Z)

res <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-NSCLC-weighted.csv",
    row.names = 1))
sivcm.parameter <- res[which.min(res[, ncol(res)]), ]
sivcm.beta <- sivcm.parameter[1:q]
sivcm.coef <- sivcm.parameter[-(1:q)]
sivcm.coef <- sivcm.coef[-length(sivcm.coef)]

nknot <- 3
degree <- 2
d <- nknot - 2 + degree
knots <- qnorm(seq(0, 1, l = nknot)[-c(1, nknot)])
Boundary.knots <- c(-max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))),
                    max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))))

B <- bSpline.ori(U %*% sivcm.beta, knots = knots, Boundary.knots = Boundary.knots)$B
sivcm.g <- get.g(X, U, Z, B, family = "gaussian", sivcm.beta, sivcm.coef)
colnames(sivcm.g) <- c("Intercept", colnames(X), colnames(Z), "B1", "B2", "B3")
sivcm.coef.mat <- matrix(sivcm.coef[(1:(p * (d + 1))) + 1], ncol = d + 1, byrow = TRUE)
rownames(sivcm.coef.mat) <- colnames(X)
sel.vary <- which(sivcm.coef.mat[, 2] != 0) + 1
for (i in sel.vary) {
    assign(paste("est.g", i, sep = ""), approxfun(U %*% sivcm.beta, sivcm.g[, i]))
}

gridpt <- seq(min(U %*% sivcm.beta), max(U %*% sivcm.beta), l = 100)
pdf("./RealDataAnalysis/RealDataAnalysisResults/Figure1.pdf")
par(mfrow = c(3, 3))
for (i in sel.vary) {
    plot(gridpt, do.call(paste("est.g", i, sep = ""), list(gridpt)),
         main = colnames(sivcm.g)[i], lty = 1, type = "l", lwd = 2, xlab = "", ylab = "",
         ylim = c(min(do.call(paste("est.g", i, sep = ""), list(gridpt))) - 5,
                  max(do.call(paste("est.g", i, sep = ""), list(gridpt))) + 5))
}
plot(gridpt, bSpline.ori(gridpt, Boundary.knots = Boundary.knots)$B %*% 
         tail(sivcm.coef, 3) + sivcm.coef[1], type = "l", main = "Intercept",
     xlab = "", ylab = "", lty = 1, lwd = 2,
     ylim = c(min(bSpline.ori(gridpt, Boundary.knots = Boundary.knots)$B %*%
                      tail(sivcm.coef, 3) + sivcm.coef[1]) - 2, 
              max(bSpline.ori(gridpt, Boundary.knots = Boundary.knots)$B %*%
                      tail(sivcm.coef, 3) + sivcm.coef[1]) + 2))
dev.off()

# Estimated index parameter
sivcm.beta

# Estimated coefficients
sivcm.coef.mat[apply(sivcm.coef.mat, 1, function(x) any(x != 0)),]


#-------------------------------------------------------------------------------
#                                    Figure 3
#-------------------------------------------------------------------------------

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

res <- as.matrix(read.csv(
    "./RealDataAnalysis/RealDataAnalysisResults/RealDataAnalysis-LGG-weighted.csv",
    row.names = 1))
sivcm.parameter <- res[which.min(res[, ncol(res)]), ]
sivcm.beta <- sivcm.parameter[1:q]
sivcm.coef <- sivcm.parameter[-(1:q)]
sivcm.coef <- sivcm.coef[-length(sivcm.coef)]

nknot <- 3
degree <- 2
d <- nknot - 2 + degree
knots <- qnorm(seq(0, 1, l = nknot)[-c(1, nknot)])
Boundary.knots <- c(-max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))), 
                    max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))))

B <- bSpline.ori(U %*% sivcm.beta, knots = knots, Boundary.knots = Boundary.knots)$B
sivcm.g <- get.g(X, U, Z, B, family = "cox", sivcm.beta, sivcm.coef)
colnames(sivcm.g) <- c(colnames(X), colnames(Z), "B1", "B2", "B3")
sivcm.coef.mat <- matrix(sivcm.coef[(1:(p * (d + 1)))], ncol = d + 1, byrow = TRUE)
rownames(sivcm.coef.mat) <- colnames(X)
sel.vary <- which(sivcm.coef.mat[, 2] != 0)
for (i in sel.vary) assign(paste("est.g", i, sep = ""), approxfun(U %*% sivcm.beta, sivcm.g[, i]))

gridpt <- seq(min(U %*% sivcm.beta), max(U %*% sivcm.beta), l = 100)
pdf("./RealDataAnalysis/RealDataAnalysisResults/Figure2.pdf")
par(mfrow = c(3, 3))
for (i in sel.vary) {
    plot(gridpt, do.call(paste("est.g", i, sep = ""), list(gridpt)),
         main = "Cyclin B1", lty = 1, type = "l", lwd = 2, xlab = "", ylab = "",
         ylim = c(min(do.call(paste("est.g", i, sep = ""), list(gridpt))) - 1, 
                  max(do.call(paste("est.g", i, sep = ""), list(gridpt))) + 1))
}
plot(gridpt, bSpline.ori(gridpt, Boundary.knots = Boundary.knots)$B %*% 
         tail(sivcm.coef, 3) + sivcm.coef[1], type = "l", main = "Intercept",
     xlab = "", ylab = "", lty = 1, lwd = 2,
     ylim = c(min(bSpline.ori(gridpt, Boundary.knots = Boundary.knots)$B %*%
                      tail(sivcm.coef, 3) + sivcm.coef[1]) - 1, 
              max(bSpline.ori(gridpt, Boundary.knots = Boundary.knots)$B %*%
                      tail(sivcm.coef, 3) + sivcm.coef[1]) + 1))
dev.off()

# Estimated index parameter
sivcm.beta

# Estimated coefficients
sivcm.coef.mat[apply(sivcm.coef.mat, 1, function(x) any(x != 0)), ]
