## FRM HW8
# Part 1
require(Rsolnp)
require(Rglpk)
# (1) estimate covariance matrices in 2 regimes
retData <- read.csv("~/R/ORF535/hw8_regime.csv", header = TRUE, stringsAsFactors = FALSE)
retData_norm <- retData[which(retData$Regime == 1), ]
retData_crash <- retData[which(retData$Regime == -1), ]

Q_norm <- matrix(0, nrow = 4, ncol = 4)
Q_norm[1:3, 1:3] <- cov(retData_norm[ , 2:4]) * 256
Q_norm[4, 4] <- 0.0001
Q_crash <- matrix(0, nrow = 4, ncol = 4)
Q_crash[1:3, 1:3] <- cov(retData_crash[ , 2:4]) * 256
Q_crash[4, 4] <- 0.0001

# (2) simulate 100 scenarios
r_norm <- c(0.04, 0.09, 0.025, 0.002)
r_crash <- c(0.07, -0.17, 0.1, 0)
expectedRet <- r_norm * 0.85 + r_crash * 0.15

simRet <- matrix(0, ncol = 4, nrow = 100)
for (i in 1:85)  simRet[i, ] <- 1 + r_norm + t(rnorm(4)) %*% chol(Q_norm)
for (i in 1:15)  simRet[(i + 85), ] <- 1 + r_crash + t(rnorm(4)) %*% chol(Q_crash)

# (3) optimization
vec <- seq(102, 107, by = 1)
res_part1 <- matrix(0, nrow = length(vec), ncol = 7)
colnames(res_part1) <- c("SP500", "Bonds", "Cash", "FTSE", "VaR", "cVaR", "Expected Wealth")

for (i in 1:length(vec)){
  opt.fun <- function(theta){
    bVaR <- (100 - (100 * simRet[, 2:4] %*% theta[1:3] + theta[4] * (simRet[, 1] - 1))) - theta[5]
    bVaR <- cbind(bVaR, rep(0, nrow(simRet)))
    return(theta[5] + 1 / ((1 - 0.95) * nrow(simRet)) * sum(apply(bVaR, 1, max)))
  }
  eq.fun <- function(theta){
    theta[1] + theta[2] + theta[3]
  }
  ineq.fun <- function(theta){
    100 * (sum(theta[1:3] * expectedRet[2:4]) + 1) + theta[4] * expectedRet[1]
  }
  res <- solnp(pars = c(1/3, 1/3, 1/3, 5, 1), fun = opt.fun, eqfun = eq.fun, eqB = 1, 
               ineqfun = ineq.fun, ineqLB = vec[i], ineqUB = 500, 
               LB = rep(0, 5), UB = c(1, 1, 1, 75, 100))
  
  res_part1[i, 1:5] <- res$pars
  res_part1[i, 6] <- res$values[length(res$values)]
  res_part1[i, 7] <- ineq.fun(res$pars)
}

plot(res_part1[, 6], res_part1[, 7], type = "l", main = "", xlab = "CVaR (million)", ylab = "Expected Wealth (million)", ylim = c(103, 107))
write.csv(res_part1, "~/R/ORF535/hw8_part1.csv")



# Part 2
smart.round <- function(x) {
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}

loss <- as.matrix(read.csv("~/R/ORF535/hw8_lossMatrix.csv", header = FALSE, stringsAsFactors = FALSE))
prob <- unlist(read.csv("~/R/ORF535/hw8_prob.csv", header = FALSE, stringsAsFactors = FALSE))
netRevenue <- unlist(read.csv("~/R/ORF535/hw8_netRevenue.csv", header = FALSE, stringsAsFactors = FALSE))
expectedLoss <- as.vector(t(prob) %*% loss)

# simulate scenarios for loss 
simLoss <- matrix(0, ncol = 7, nrow = 1000)
numScen <- smart.round(1000 * prob)
rowCount <- 0

for (i in 1:length(prob)){
  for (j in 1:numScen[i]) simLoss[(rowCount + j), ] <- loss[i, ]
  rowCount <- rowCount + numScen[i]
}
scen <- cbind(simRet[rep(c(1:100), rep(1000, 100)), ], simLoss[rep(c(1:1000), 100), ])


# (2) optimization
surplus <- function(theta){
  capital0 <- 100 + theta[8] + sum(theta[1:7] * netRevenue)
  r <- scen[, 2:4] %*% theta[9:11]
  l <- scen[, 5:11] %*% theta[1:7]
  surplus <- capital0 * r + theta[12] * (scen[, 1] - 1) - l - 1.04 * theta[8]
  return(surplus)
}

opt.fun2 <- function(theta){
  capital0 <- 100 + theta[8] + sum(theta[1:7] * netRevenue)
  r <- scen[, 2:4] %*% theta[9:11]
  l <- scen[, 5:11] %*% theta[1:7]
  surplus <- capital0 * r + theta[12] * (scen[, 1] - 1) - l - 1.04 * theta[8]
  
  bVaR2 <- (100 - as.vector(surplus)) - theta[13]
  bVaR2 <- cbind(bVaR2, rep(0, length(bVaR2)))
  return(theta[13] + 1 / ((1 - 0.99) * length(bVaR2)) * sum(apply(bVaR2, 1, max)))
}
eq.fun2 <- function(theta){
  theta[9] + theta[10] + theta[11]
}
ineq.fun2 <- function(theta){
  (100 + theta[8] +  sum(theta[1:7] * netRevenue)) * (sum(theta[9:11] * expectedRet[2:4]) + 1) + theta[12] * expectedRet[1] - sum(theta[1:7] * expectedLoss) - 1.04 * theta[8]
}

vec_part2 <- seq(105, 125, by = 5)
res_part2 <- matrix(0, nrow = length(vec_part2), ncol = 18)
colnames(res_part2) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", 
                         "Borrowing", "SP500", "Bonds", "Cash", "FTSE", 
                         "Expected Surplus", "VaR", "cVaR", "Prob of Default", 
                         "Expected Deficit", "AssetToCapital")

for (i in 1:length(vec_part2)){
  res2 <- solnp(pars = c(rep(0.5, 7), 20, rep(1/3, 3), 20, 1), fun = opt.fun2, eqfun = eq.fun2, eqB = 1, 
                ineqfun = ineq.fun2, ineqLB = vec_part2[i], ineqUB = 500, 
                LB = rep(0, 13), UB = c(rep(1, 7), 50, rep(1, 3), 75, 500))
  
  res_part2[i, 1:12] <- res2$pars[1:12]
  res_part2[i, 13] <- ineq.fun2(res2$pars)
  res_part2[i, 14] <- res2$pars[13]
  res_part2[i, 15] <- res2$values[length(res2$values)]
  
  S <- surplus(res2$pars)
  res_part2[i, 16] <- sum(S < 0)/length(S)
  res_part2[i, 17] <- mean(S[S < 0])
  
  A <- 100 + res2$pars[8] + sum(res2$pars[1:7] * netRevenue)
  C <- 100 + sum(res2$pars[1:7] * expectedLoss) * 0.15
  res_part2[i, 18] <- A / C
}



plot(res_part2[, 15], res_part2[, 13], type = "l", main = "", xlab = "1%-CVaR (million)", ylab = "Expected Surplus (million)", ylim = c(100, 130))
write.csv(res_part2, "~/R/ORF535/hw8_part2.csv")

# Part 3
opt.fun3 <- function(theta){
  
  capital0 <- 100 + theta[8] + sum(theta[1:7] * netRevenue)
  r <- scen[, 2:4] %*% theta[9:11]
  l <- scen[, 5:11] %*% theta[1:7]
  surplus <- capital0 * r + theta[12] * (scen[, 1] - 1) - l - 1.04 * theta[8]
  
  
  bVaR2 <- (100 - as.vector(surplus)) - theta[13]
  bVaR2 <- cbind(bVaR2, rep(0, length(bVaR2)))
  return(theta[13] + 1 / ((1 - 0.99) * length(bVaR2)) * sum(apply(bVaR2, 1, max)) + 1000 * sum(theta[1:7]*(1 - theta[1:7])))
}
eq.fun3 <- function(theta){
  theta[9] + theta[10] + theta[11]
}
ineq.fun3 <- function(theta){
  (100 + theta[8] +  sum(theta[1:7] * netRevenue)) * (sum(theta[9:11] * expectedRet[2:4]) + 1) + theta[12] * expectedRet[1] - sum(theta[1:7] * expectedLoss) - 1.04 * theta[8]
}


vec_part3 <- seq(105, 125, by = 5)
res_part3 <- matrix(0, nrow = length(vec_part3), ncol = 18)
colnames(res_part3) <- c("y1", "y2", "y3", "y4", "y5", "y6", "y7", 
                         "Borrowing", "SP500", "Bonds", "Cash", "FTSE", 
                         "Expected Surplus", "VaR", "cVaR", "Prob of Default", 
                         "Expected Deficit", "AssetToCapital")


for (i in 1:length(vec_part2)){
  res3 <- solnp(pars = c(rep(0.5, 7), 20, rep(1/3, 3), 20, 1), fun = opt.fun3, eqfun = eq.fun3, eqB = 1, 
                ineqfun = ineq.fun3, ineqLB = vec_part3[i], ineqUB = 500, 
                LB = rep(0, 13), UB = c(rep(1, 7), 50, rep(1, 3), 75, 500))
  
  res_part3[i, 1:12] <- res3$pars[1:12]
  res_part3[i, 13] <- ineq.fun3(res3$pars)
  res_part3[i, 14] <- res3$pars[13]
  res_part3[i, 15] <- res3$values[length(res3$values)]
  
  S <- surplus(res3$pars)
  res_part3[i, 16] <- sum(S < 0)/length(S)
  res_part3[i, 17] <- mean(S[S < 0])
  
  A <- 100 + res3$pars[8] + sum(res3$pars[1:7] * netRevenue)
  C <- 100 + sum(res3$pars[1:7] * expectedLoss) * 0.15
  res_part3[i, 18] <- A / C
}

plot(res_part3[, 15], res_part3[, 13], type = "l",main = "", xlab = "1%-CVaR (million)", ylab = "Expected Surplus (million)", ylim = c(100, 130))
write.csv(res_part3, "~/R/ORF535/hw8_part3.csv")
