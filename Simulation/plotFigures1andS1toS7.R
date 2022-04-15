#-------------------------------------------------------------------------------
#                    Plot figures for the simulation studies
#-------------------------------------------------------------------------------

# Figure 1: family = "gaussian" and p = 100
# Figure S1: family = "gaussian" and p = 20
# Figure S2: family = "gaussian" and p = 50
# Figure S3: family = "gaussian" and p = 300
# Figure S4: family = "cox" and p = 20
# Figure S5: family = "cox" and p = 50
# Figure S6: family = "cox" and p = 100
# Figure S7: family = "cox" and p = 300

library(psivcm)

model <- 1

for(family in c("gaussian", "cox")){
  for(p in c(20, 50, 100, 300)){
    for (initbeta in 1:5) {
      assign(paste("res", initbeta, ".unweighted", sep = ""), 
             as.matrix(read.csv(paste(
               "./Simulation/SimulationResults/SimulationResults-setting",
               model, "-", family, "-p", p, "-beta", initbeta, 
               "-unweighted.csv", sep = ""), row.names = 1)))   
    }
    for (initbeta in 1:5) {
      assign(paste("est", initbeta, ".unweighted", sep = ""), 
             as.matrix(read.csv(paste(
               "./Simulation/SimulationResults/Estimates-setting",
               model, "-", family, "-p", p, "-beta", initbeta,
               "-unweighted.csv", sep = ""))))   
    }
    
    res1.unweighted[(res1.unweighted[,"converge"] == 0), "ic"] <- 1e+10
    res2.unweighted[(res2.unweighted[,"converge"] == 0), "ic"] <- 1e+10
    res3.unweighted[(res3.unweighted[,"converge"] == 0), "ic"] <- 1e+10
    res4.unweighted[(res4.unweighted[,"converge"] == 0), "ic"] <- 1e+10
    res5.unweighted[(res5.unweighted[,"converge"] == 0), "ic"] <- 1e+10
    
    initbeta.opt.unweighted <- apply(cbind(res1.unweighted[, "ic"], 
                                           res2.unweighted[, "ic"], 
                                           res3.unweighted[, "ic"], 
                                           res4.unweighted[, "ic"], 
                                           res5.unweighted[, "ic"]), 
                                     1, which.min)
    
    for (initbeta in 1:5) {
      assign(paste("res", initbeta, ".weighted", sep = ""), 
             as.matrix(read.csv(paste(
               "./Simulation/SimulationResults/SimulationResults-setting",
               model, "-", family, "-p", p, "-beta", initbeta,
               "-weighted.csv", sep = ""), row.names = 1)))   
    }
    for (initbeta in 1:5) {
      assign(paste("est", initbeta, ".weighted", sep = ""), 
             as.matrix(read.csv(paste(
               "./Simulation/SimulationResults/Estimates-setting",
               model, "-", family, "-p", p, "-beta", initbeta,
               "-weighted.csv", sep = ""))))
    }
    
    res1.weighted[(res1.weighted[,"converge"] == 0), "ic"] <- 1e+10
    res2.weighted[(res2.weighted[,"converge"] == 0), "ic"] <- 1e+10
    res3.weighted[(res3.weighted[,"converge"] == 0), "ic"] <- 1e+10
    res4.weighted[(res4.weighted[,"converge"] == 0), "ic"] <- 1e+10
    res5.weighted[(res5.weighted[,"converge"] == 0), "ic"] <- 1e+10
    
    initbeta.opt.weighted <- apply(cbind(res1.weighted[, "ic"], 
                                         res2.weighted[, "ic"],
                                         res3.weighted[, "ic"], 
                                         res4.weighted[, "ic"], 
                                         res5.weighted[, "ic"]), 
                                   1, which.min)
    
    result.sivcm.unweighted <- matrix(nrow = 100, ncol = ncol(res1.unweighted))
    for (seed in 1:100) {
      if (initbeta.opt.unweighted[seed] == 1) {
        result.sivcm.unweighted[seed, ] <- res1.unweighted[seed, ]
      } else if (initbeta.opt.unweighted[seed] == 2) {
        result.sivcm.unweighted[seed, ] <- res2.unweighted[seed, ]    
      } else if (initbeta.opt.unweighted[seed] == 3) {
        result.sivcm.unweighted[seed, ] <- res3.unweighted[seed, ]    
      } else if (initbeta.opt.unweighted[seed] == 4) {
        result.sivcm.unweighted[seed, ] <- res4.unweighted[seed, ]
      } else if (initbeta.opt.unweighted[seed] == 5) {
        result.sivcm.unweighted[seed, ] <- res5.unweighted[seed, ]        
      }   
    }
    result.sivcm.weighted <- matrix(nrow = 100, ncol = ncol(res1.weighted))
    for (seed in 1:100) {
      if (initbeta.opt.weighted[seed] == 1) {
        result.sivcm.weighted[seed, ] <- res1.weighted[seed, ]
      } else if (initbeta.opt.weighted[seed] == 2) {
        result.sivcm.weighted[seed, ] <- res2.weighted[seed, ]
      } else if (initbeta.opt.weighted[seed] == 3) {
        result.sivcm.weighted[seed, ] <- res3.weighted[seed, ]
      } else if (initbeta.opt.weighted[seed] == 4) {
        result.sivcm.weighted[seed, ] <- res4.weighted[seed, ]
      } else if (initbeta.opt.weighted[seed] == 5) {
        result.sivcm.weighted[seed, ] <- res5.weighted[seed, ]
      }
    }
    
    store.xlim <- matrix(nrow = 100, ncol = 2)
    for (seed in 1:100) {
      
      data.train <- as.matrix(read.csv(paste(
        "./Simulation/SimulationData/SimulationData-p", p, "-", seed, 
        ".csv", sep = "")))
      if (family == "cox")
        Y <- Surv(data.train[, "time"], data.train[, "status"])
      if (family == "gaussian")
        Y <- data.train[, "Y"]
      X <- data.train[, which(sapply(1:ncol(data.train), function(x) 
        unlist(strsplit(colnames(data.train)[x], ""))[1]) == "X")]
      U <- data.train[, which(sapply(1:ncol(data.train), function(x) 
        unlist(strsplit(colnames(data.train)[x], ""))[1]) == "U")]
      Z <- U[, -4]
      
      if (initbeta.opt.unweighted[seed] == 1) {
        sivcm.est.unweighted <- est1.unweighted[seed, ]
      } else if (initbeta.opt.unweighted[seed] == 2) {
        sivcm.est.unweighted <- est2.unweighted[seed, ]
      } else if (initbeta.opt.unweighted[seed] == 3) {
        sivcm.est.unweighted <- est3.unweighted[seed, ]
      } else if (initbeta.opt.unweighted[seed] == 4) {
        sivcm.est.unweighted <- est4.unweighted[seed, ]
      } else if (initbeta.opt.unweighted[seed] == 5) {
        sivcm.est.unweighted <- est5.unweighted[seed, ]
      }
      sivcm.beta.unweighted <- sivcm.est.unweighted[1:4]
      sivcm.coef.unweighted <- sivcm.est.unweighted[- (1:4)]
      
      if (initbeta.opt.weighted[seed] == 1) {
        sivcm.est.weighted <- est1.weighted[seed, ]
      } else if (initbeta.opt.weighted[seed] == 2) {
        sivcm.est.weighted <- est2.weighted[seed, ]
      } else if (initbeta.opt.weighted[seed] == 3) {
        sivcm.est.weighted <- est3.weighted[seed, ]
      } else if (initbeta.opt.weighted[seed] == 4) {
        sivcm.est.weighted <- est4.weighted[seed, ]
      } else if (initbeta.opt.weighted[seed] == 5) {
        sivcm.est.weighted <- est5.weighted[seed, ]
      }
      sivcm.beta.weighted <- sivcm.est.weighted[1:4]
      sivcm.coef.weighted <- sivcm.est.weighted[-(1:4)]
      
      
      store.xlim[seed, 1] <- max(min(U %*% sivcm.beta.unweighted), 
                                 min(U %*% sivcm.beta.weighted))
      store.xlim[seed, 2] <- min(max(U %*% sivcm.beta.unweighted), 
                                 max(U %*% sivcm.beta.weighted))
      
      knots <- qnorm(seq(0, 1, l = 3)[- c(1, 3)])
      Boundary.knots <- c(-max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))), 
                          max(apply(U, 1, function(x) sqrt(sum(x ^ 2)))))
      
      if (sivcm.beta.unweighted[1] + sivcm.beta.unweighted[3] + 
          sivcm.beta.unweighted[4] - sivcm.beta.unweighted[2] < 0) {
        B <- bSpline.ori(U %*% -sivcm.beta.unweighted, knots = knots,
                         Boundary.knots = Boundary.knots)$B
      } else {
        B <- bSpline.ori(U %*% sivcm.beta.unweighted, knots = knots,
                         Boundary.knots = Boundary.knots)$B
      }
      
      if (family == "cox") {
        for (i in 1:20) {
          assign(paste("est.g.unweighted", i, ".seed", seed, sep = ""),
                 approxfun(U %*% sivcm.beta.unweighted, 
                           get.g(X, U, Z, B, sivcm.beta.unweighted,
                                 sivcm.coef.unweighted, family = family)[, i]))
        }
      }
      if (family == "gaussian") {
        for (i in 1:20) {
          assign(paste("est.g.unweighted", i, ".seed", seed, sep = ""),
                 approxfun(U %*% sivcm.beta.unweighted, 
                           get.g(X, U, Z, B, sivcm.beta.unweighted,
                                 sivcm.coef.unweighted, family = family)[, i + 1]))
        }
      }
      if (sivcm.beta.weighted[1] + sivcm.beta.weighted[3] + 
          sivcm.beta.weighted[4] - sivcm.beta.weighted[2] < 0) {
        B <- bSpline.ori(U %*% -sivcm.beta.weighted, knots = knots, 
                         Boundary.knots = Boundary.knots)$B
      } else {
        B <- bSpline.ori(U %*% sivcm.beta.weighted, knots = knots, 
                         Boundary.knots = Boundary.knots)$B
      }
      
      if (family == "cox") {
        for (i in 1:20) {
          assign(paste("est.g.weighted", i, ".seed", seed, sep = ""),
                 approxfun(U %*% sivcm.beta.weighted, 
                           get.g(X, U, Z, B, sivcm.beta.weighted, 
                                 sivcm.coef.weighted, family = family)[, i]))
        }
      }
      if (family == "gaussian") {
        for (i in 1:20) {
          assign(paste("est.g.weighted", i, ".seed", seed, sep = ""),
                 approxfun(U %*% sivcm.beta.weighted, 
                           get.g(X, U, Z, B, sivcm.beta.weighted, 
                                 sivcm.coef.weighted, family = family)[, i + 1]))
        }
        
      }
    }
    
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
        g19 <- function(z) 1 / (1 + exp(-z * 5))
        g20 <- function(z) 1 / (1 + exp(z * 4)) - 0.3
      }
    }
    
    average.grid.unweighted <- matrix(nrow = 20, ncol = 100)
    average.grid.weighted <- matrix(nrow = 20, ncol = 100)
    gridpt <- seq(max(store.xlim[, 1]), min(store.xlim[, 2]), l = 100)
    for (i in 1:20) {
      sum <- 0
      for (seed in 1:100) {
        sum <- sum + do.call(paste("est.g.unweighted", i, ".seed", seed, sep = ""),
                             list(gridpt))
      }
      average.grid.unweighted[i, ] <- sum / 100
    }
    for (i in 1:20) {
      sum <- 0
      for (seed in 1:100) {
        sum <- sum + do.call(paste("est.g.weighted", i, ".seed", seed, sep = ""),
                             list(gridpt))
      }
      average.grid.weighted[i, ] <- sum / 100
    }
    
    pdf(paste("./Simulation/SimulationResults/SimulationResults-setting1-", 
              family, "-p", p, "-coef.pdf", sep = ""))
    par(oma = c(2, 2, 2, 2), mar = c(2, 2, 2, 2))
    
    m <- matrix(c(1:20, 21, 21, 21, 21), nrow = 6, ncol = 4, byrow = TRUE)
    layout(mat = m, heights = c(0.14, 0.14, 0.14, 0.14, 0.14, 0.1))
    
    for (i in 1:20) {
      xlim <- c(range(gridpt)[1] * 1.1, range(gridpt)[2] * 1.1)
      ylim <- c(-2.4, 2.4)
      plot(gridpt, do.call(paste("g", i, sep = ""), list(gridpt)), 
           ylim = ylim, xlim = xlim, main = bquote(italic(g)[.(i)]), 
           lty = 1, type = "l",
           ylab = paste("g", i, "(grid point)", sep = ""),
           xlab = "grid point", lwd = 1, cex.axis = 0.9)
      lines(gridpt, average.grid.unweighted[i, ], col = "green", lwd = 2, lty = 2)
      lines(gridpt, average.grid.weighted[i, ], col = "red", lwd = 2, lty = 3)
    }
    
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("top", legend = c("True", "Unweighted", "Weighted"), 
           col = c("black", "green", "red"), horiz = TRUE, lwd = c(1, 2, 2), lty = 1:3)
    
    dev.off()
  }
  
}
