library(SISIR)

## PACKAGE
library(DiceKriging)

################################################################################
## Generation exemple
set.seed(118)
n <- 100
ntest <- 30
H <- 10
p <- 200
tmax <- 1.3
tslice <- c(.2, .4, .8, 1.1)
pb <- 1
relevant.slices <- list(1)

################################################################################
# Generate data
source("generate_data.R")
op <- par()
dat <- generate_data(pb = pb, ntrain = n, ntest = ntest, p = p, draw=TRUE, 
                     tslice = tslice, relevant.slices = relevant.slices)
X      <- dat$X
Y      <- dat$Y

################################################################################
#--- Ridge estimation ---------------------------------------
## Selection of d and mu_2
list_mu2 <- 10^(-3:5)
listH <- c(10)
list_d <- 1:4
res_tune <- tune.ridgeSIR(X, Y, listH, list_mu2, list_d, nfolds = 10,
                          parallel = FALSE)

printable_res <- matrix(res_tune[ ,"cverror"], ncol = length(list_d),
                        byrow = TRUE)
rownames(printable_res) <- log10(list_mu2)
colnames(printable_res) <- list_d
par(op)
plot(log10(list_mu2), printable_res[ ,ncol(printable_res)], main = "",
     ylab = "CV error", xlab = expression(log[10](mu[2])), type = "b",
     pch = "+")
# dev.print(png, file = "../results/cverror-test1.png", width = 600, height = 500, res = 150)
chosen_mu <- which.min(printable_res[ ,ncol(printable_res)])
### chosen_mu: 1

### Selection of d
printable_res2 <- matrix(res_tune[res_tune$H == 10, "Rd"], 
                         ncol = length(list_d), byrow = TRUE)
rownames(printable_res2) <- log10(list_mu2)
colnames(printable_res2) <- list_d
par(op)
plot(list_d, printable_res2[chosen_mu, ], type = "b", pch = "+", xlab = "d",
     ylab = expression(R(d)))
# dev.print(png, file = "../results/Rd-test1.png", width = 600, height = 500, res = 150)
d <- 1

### Final check
par(op)
plot(log10(list_mu2), printable_res[ ,d], main = "", ylab = "CV error",
     xlab = expression(log[10](mu[2])), type = "b", pch = "+")
which.min(printable_res[ ,ncol(printable_res)]) ## consistent
                  
### Ridge step  
res_ridge <- ridgeSIR(X, Y, H = 10, d = 1, mu2 = list_mu2[chosen_mu])

################################################################################
#--- Sparse estimation --------------------------------------
res_sparse <- SISIR(res_ridge, parallel = FALSE)

## analysis of results
bestBIC <- which.min(res_sparse$quality$BIC)
bestAIC <- which.min(res_sparse$quality$AIC)
bestCV <- which.min(res_sparse$quality$CV)

par(mfrow=c(1,2))
plot(res_sparse$quality$nbint, res_sparse$quality$BIC, type = "l",
     ylim = range(c(res_sparse$quality$BIC, res_sparse$quality$AIC)))
lines(res_sparse$quality$nbint, res_sparse$quality$AIC, col = "red")
legend("topleft", col=c("red","black"), legend=c("AIC","BIC"), lty=c(1,1,1))

plot(res_sparse$quality$nbint, res_sparse$quality$CVerror, col="red", type="l")

## plots
grid <- seq(0.2, 1.1, length=p)
active_interval <- ((grid >= 0.2) & (grid <= 0.4))
r1 <- which((tslice >= 0.2) & (tslice <= 0.4))

### best model
par(op)
plot(grid, res_sparse$sEDR[[bestCV]], ylab = expression(paste(hat(a))[1]),
     type = "l", xlab = "")
lines(grid[active_interval], res_sparse$sEDR[[bestCV]][active_interval],
      col = "red")
# dev.print(png, file = "../results/beta_M1_best.png", width = 600, height = 500, res = 150)

### p=200 intervals (usual LASSO estimate)
plot(grid, res_sparse$sEDR[[1]], ylab = expression(paste(hat(a))[1]),
     type = "l", xlab = "")
lines(grid[active_interval], res_sparse$sEDR[[1]][active_interval], col="red")
# dev.print(png, file="../results/beta_M1_200int.png", width=600, height=500, res=150)

### 41 intervals
plot(grid, res_sparse$sEDR[[38]], ylab = expression(paste(hat(a))[1]),
     type = "l",xlab = "")
lines(grid[active_interval], res_sparse$sEDR[[38]][active_interval], col="red")
# dev.print(png, file = "../results/beta_M1_41int.png", width = 600, height = 500, res = 150)

## 5 intervals
plot(grid, res_sparse$sEDR[[79]], ylab=expression(paste(hat(a))[1]),
     type = "l", xlab = "")
lines(grid[active_interval], res_sparse$sEDR[[79]][active_interval], col="red")
# dev.print(png, file = "../results/beta_M1_5int.png", width = 600, height = 500, res = 150)
