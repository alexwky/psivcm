#-------------------------------------------------------------------------------
#                  Generate tables for the simulation results
#-------------------------------------------------------------------------------

output <- list()
output[[1]] <- c("sen", "sen.int", "fdr", "fdr.int", "card", "card.int", "aip",
                 "mse", "c.index")
output[[2]] <- c("sen", "fdr", "card", "card.main", "card.int", "mse",
                 "c.index")
output[[3]] <- c("sen", "sen.int", "fdr", "fdr.int", "card", "card.int", "aip",
                 "mse", "c.index")
output[[4]] <- c("sen", "sen.int", "fdr", "fdr.int", "card", "card.int", "aip",
                 "mse", "c.index")

extract <- function(model, family, tune, nonlinear.main, p.seq){
  
  sivc <- matrix(nrow = length(output[[model]]), ncol = length(p.seq))
  colnames(sivc) <- sapply(1:length(p.seq), function(x) paste("p=", p.seq[x], sep = ""))
  rownames(sivc) <- output[[model]]
  sivc.w <- matrix(nrow = length(output[[model]]), ncol = length(p.seq))
  colnames(sivc.w) <- sapply(1:length(p.seq), function(x) paste("p=", p.seq[x], sep = ""))
  rownames(sivc.w) <- output[[model]]
  main <- matrix(nrow = length(output[[model]]), ncol = length(p.seq))
  colnames(main) <- sapply(1:length(p.seq), function(x) paste("p=", p.seq[x], sep = ""))
  rownames(main) <- output[[model]]
  main.w <- matrix(nrow = length(output[[model]]), ncol = length(p.seq))
  colnames(main.w) <- sapply(1:length(p.seq), function(x) paste("p=", p.seq[x], sep = ""))
  rownames(main.w) <- output[[model]]
  int <- matrix(nrow = length(output[[model]]), ncol = length(p.seq))
  colnames(int) <- sapply(1:length(p.seq), function(x) paste("p=", p.seq[x], sep = ""))
  rownames(int) <- output[[model]]
  int.w <- matrix(nrow = length(output[[model]]), ncol = length(p.seq))
  colnames(int.w) <- sapply(1:length(p.seq), function(x) paste("p=", p.seq[x], sep = ""))
  rownames(int.w) <- output[[model]]
  
  for(p in p.seq){
    for (initbeta in 1:5) {
      assign(paste("res", initbeta, ".unweighted", sep = ""), 
             as.matrix(read.csv(paste(
               "./Simulation/SimulationResults/SimulationResults-setting",
               model, ifelse(nonlinear.main == TRUE, "-", 
                             ifelse(model == 3, "-mis-", "-")),
               family, "-p", p, "-beta", initbeta, "-unweighted",
               ifelse(tune == "aic", "-aic.csv", ".csv"), sep = ""),
               row.names = 1)))   
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
                                           res5.unweighted[, "ic"]), 1, which.min)
    for (initbeta in 1:5) {
      assign(paste("res", initbeta, ".weighted", sep = ""), 
             as.matrix(read.csv(paste(
               "./Simulation/SimulationResults/SimulationResults-setting",
               model, ifelse(nonlinear.main == TRUE, "-", 
                             ifelse(model == 3, "-mis-", "-")),
               family, "-p", p, "-beta", initbeta, "-weighted",
               ifelse(tune == "aic", "-aic.csv", ".csv"), sep = ""),
               row.names = 1)))   
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
                                         res5.weighted[, "ic"]), 1, which.min)
    
    result.sivcm.unweighted <- matrix(nrow = 100, ncol = ncol(res1.unweighted))
    colnames(result.sivcm.unweighted) <- colnames(res1.unweighted)
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
    colnames(result.sivcm.weighted) <- colnames(res1.weighted)
    
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
    sivc[, which(p == p.seq)] <- colMeans(result.sivcm.unweighted)[output[[model]]]
    sivc.w[, which(p == p.seq)] <- colMeans(result.sivcm.weighted)[output[[model]]]
    
    if(model != 3){
      main[, which(p == p.seq)] <- colMeans(as.matrix(read.csv(
        paste("./Simulation/SimulationResults/SimulationResults-setting",
              model, "-", family, "-p", p, "-main-unweighted", 
              ifelse(tune == "aic", "-aic.csv", ".csv"),
              sep = "")))[, -1])[output[[model]]]
      main.w[, which(p == p.seq)] <- colMeans(as.matrix(read.csv(
        paste("./Simulation/SimulationResults/SimulationResults-setting",
              model, "-", family, "-p", p, "-main-weighted",
              ifelse(tune == "aic", "-aic.csv", ".csv"),
              sep = "")))[, -1])[output[[model]]]
      int[, which(p == p.seq)] <- colMeans(as.matrix(read.csv(
        paste("./Simulation/SimulationResults/SimulationResults-setting",
              model, "-", family, "-p", p, "-int-unweighted",
              ifelse(tune == "aic", "-aic.csv", ".csv"),
              sep = "")))[, -1])[output[[model]]]
      int.w[, which(p == p.seq)] <- colMeans(as.matrix(read.csv(
        paste("./Simulation/SimulationResults/SimulationResults-setting",
              model, "-", family, "-p", p, "-int-weighted",
              ifelse(tune == "aic", "-aic.csv", ".csv"),
              sep = "")))[, -1])[output[[model]]]
    }
  }
  
  list(SIVC = sivc, SIVC.w = sivc.w, MAIN = main, MAIN.w = main.w, INT = int, INT.w = int.w)
}


# Table 1:
extract(model = 1, family = "gaussian", tune = "bic", nonlinear.main = FALSE,
        p.seq = c(20, 50, 100, 300))

# Table 2:
extract(model = 1, family = "cox", tune = "bic", nonlinear.main = FALSE,
        p.seq = c(20, 50, 100, 300))

# Table S3:
extract(model = 2, family = "gaussian", tune = "bic", nonlinear.main = FALSE,
        p.seq = c(20, 50, 100, 300))

# Table S4:
extract(model = 2, family = "cox", tune = "bic", nonlinear.main = FALSE,
        p.seq = c(20, 50, 100, 300))

# Table S5:
extract(model = 3, family = "gaussian", tune = "bic", nonlinear.main = FALSE,
        p.seq = c(20, 50, 100, 300))
extract(model = 3, family = "gaussian", tune = "bic", nonlinear.main = TRUE,
        p.seq = c(20, 50, 100, 300))

# Table S6:
extract(model = 3, family = "cox", tune = "bic", nonlinear.main = FALSE,
        p.seq = c(20, 50, 100, 300))
extract(model = 3, family = "cox", tune = "bic", nonlinear.main = TRUE,
        p.seq = c(20, 50, 100, 300))

# Table S7:
extract(model = 4, family = "gaussian", tune = "bic", nonlinear.main = FALSE,
        p.seq = c(20, 50, 100, 300))

# Results for the AIC
extract(model = 1, family = "gaussian", tune = "aic", nonlinear.main = FALSE,
        p.seq = c(20, 50, 100, 300))
extract(model = 1, family = "cox", tune = "aic", nonlinear.main = FALSE,
        p.seq = c(20, 50, 100, 300))
