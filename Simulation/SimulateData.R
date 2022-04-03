library(psivcm)
library(survival)

# Generate multiple initial values of beta for the simulation studies
initbeta <- matrix(nrow = 5, ncol = 4)
for (seed in 1:5) {
    set.seed(seed)
    beta0 <- rnorm(4)
    beta0 <- beta0 / sqrt(sum(beta0 ^ 2))
    initbeta[seed, ] <- beta0
}
write.table(initbeta, file = "./Simulation/SimulationData/Simulation-beta0.csv", 
            row.names = FALSE, col.names = FALSE, sep = ",")


# Generate random data for the simulation studies
for (p in c(20, 50, 100, 300)) {
    for (seed in 0:100) {
        data.cox <- gen.data(seed = seed, n = ifelse(seed == 0, 5000, 500), 
                             p = p, q = 4, beta = c(0.4, -0.4, 0.2, 0.8), 
                             psi = c(0.2, -0.2, 0.5, -0.5), family = "cox", 
                             exp.rate = ifelse(seed == 0, 1e-10, 0.25), model = 1)
        Y.cox <- data.cox$Y
        X.cox <- data.cox$X
        U.cox <- data.cox$U
        data.gaussian <- gen.data(seed = seed, n = ifelse(seed == 0, 5000, 500),
                                  p = p, q = 4, beta = c(0.4, -0.4, 0.2, 0.8), 
                                  psi = c(0.2, -0.2, 0.5, -0.5), 
                                  family = "gaussian", model = 1)
        Y.gaussian <- data.gaussian$Y
        X.gaussian <- data.gaussian$X
        U.gaussian <- data.gaussian$U
        if (!all(X.cox == X.gaussian) | !all(U.cox == U.gaussian)) stop()
        data <- cbind(Y.gaussian, Y.cox, X.gaussian, U.gaussian)
        colnames(data) <- c("Y", "time", "status", 
                            sapply(1:p, function(x) paste("X", x, sep = "")), 
                            sapply(1:4, function(x) paste("U", x, sep = "")))
        write.table(data, 
                    file = paste("./Simulation/SimulationData/SimulationData-", 
                                 "p", p, "-", seed, ".csv", sep = ""),
                    row.names = FALSE, col.names = TRUE, sep = ",")
    }
}


# Generate random data for the additional simulation studies
for (model in c(2, 3)) {
    for (p in c(20, 50, 100, 300)) {
        for (seed in 0:100) {
            data.cox <- gen.data(seed = seed, n = ifelse(seed == 0, 5000, 500), 
                                 p = p, q = 4, beta = c(0.4, -0.4, 0.2, 0.8), 
                                 psi = c(0.2, -0.2, 0.5, -0.5), family = "cox",
                                 exp.rate = ifelse(seed == 0, 1e-10,
                                                   ifelse(model == 2, 0.25, 0.22)), 
                                 model = model)
            Y.cox <- data.cox$Y
            X.cox <- data.cox$X
            U.cox <- data.cox$U
            data.gaussian <- gen.data(seed = seed, n = ifelse(seed == 0, 5000, 500),
                                      p = p, q = 4, beta = c(0.4, -0.4, 0.2, 0.8),
                                      psi = c(0.2, -0.2, 0.5, -0.5), 
                                      family = "gaussian", model = model)
            Y.gaussian <- data.gaussian$Y
            X.gaussian <- data.gaussian$X
            U.gaussian <- data.gaussian$U
            if (!all(X.cox == X.gaussian) | !all(U.cox == U.gaussian)) stop()
            data <- cbind(Y.gaussian, Y.cox)
            colnames(data) <- c("Y", "time", "status")
            write.table(data, file = paste(
                "./Simulation/SimulationData/AdditionalSimulationData-setting",
                model, "-p", p, "-", seed, ".csv", sep = ""), 
                row.names = FALSE, col.names = TRUE, sep = ",")
        }
    }
}

for (model in 4) {
    for (p in c(20, 50, 100, 300)) {
        for (seed in 0:100) {
            data.gaussian <- gen.data(seed = seed, n = ifelse(seed == 0, 5000, 500),
                                      p = p, q = 4, beta = c(0.4, -0.4, 0.2, 0.8),
                                      psi = c(0.2, -0.2, 0.5, -0.5), 
                                      family = "gaussian", model = model)
            data <- as.matrix(data.gaussian$Y)
            colnames(data) <- c("Y")
            write.table(data, file = paste(
                "./Simulation/SimulationData/AdditionalSimulationData-setting",
                model, "-p", p, "-", seed, ".csv", sep = ""), 
                row.names = FALSE, col.names = TRUE, sep = ",")
        }
    }
}
