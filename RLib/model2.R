library(SISIR)

## PACKAGE
library(DiceKriging)

################################################################################
## Generation exemple
set.seed(118)
n <- 100
ntest <- 30
H <- 10
p <- 300
tmax <- 1
tslice <- c(0,.1, .2, .5, .65, .78,.9, 1)
pb <- 1
relevant.slices <- list(1,4,5) 

################################################################################
# Generate data
op <- par()
source("generate_data.R")
dat <- generate_data(pb = pb, ntrain = n, ntest = ntest, p = p, draw = TRUE, 
                     tslice = tslice, relevant.slices = relevant.slices)
X      <- dat$X
Y      <- dat$Y

################################################################################
#--- Ridge estimation ---------------------------------------
## Estimation of d and mu_2
list_mu2 <- 10^(-3:5)
listH <- c(10)
list_d <- 1:4
res_tune <- tune.ridgeSIR(X, Y, listH, list_mu2, list_d, nfolds = 10)

printable_res <- matrix(res_tune[ ,"cverror"], ncol = length(list_d), 
                        byrow = TRUE)
rownames(printable_res) <- log10(list_mu2)
colnames(printable_res) <- list_d
par(op)
plot(log10(list_mu2), printable_res[ ,ncol(printable_res)], main = "",
     ylab = "CV error", xlab = expression(log[10](mu[2])), type = "b",
     pch = "+")
chosen_mu <- which.min(printable_res[ ,ncol(printable_res)])
### chosen_mu: 10

### estimation of d
printable_res2 <- matrix(res_tune[ ,"Rd"], ncol=length(list_d), byrow = TRUE)
rownames(printable_res2) <- log10(list_mu2)
colnames(printable_res2) <- list_d
plot(list_d, printable_res2[chosen_mu, ], type = "b", pch = "+", xlab = "d",
     ylab = expression(R(d)))
d <- 1

### check
plot(log10(list_mu2), printable_res[ ,d], main = "", ylab = "CV error",
     xlab = expression(log[10](mu[2])), type = "b", pch = "+")
which.min(printable_res[ ,ncol(printable_res)]) ## not consistent
chosen_mu <- which.min(printable_res[ ,d])
# dev.print(png, file="../results/cverror-test2-2.png", width=600, height=500, res=150)
### chosen_mu: 10

### re-estimation of d
plot(list_d, printable_res2[chosen_mu, ], type = "b", pch = "+", xlab = "d",
     ylab = expression(R(d)))
# dev.print(png, file = "../results/Rd-test2-2.png", width = 600, height = 500, res = 150)
d <- 1 ## consistent
                  
## Ridge step  
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
par(op)
par(mar = c(4,2,2,2))
grid <- seq(0, 1.1, length = p)
r1 <- which((grid>=0) & (grid<=0.1))
r2 <- which((grid>=0.5) & (grid<=0.65))
r3 <- which((grid>=0.65) & (grid<=0.78))
int <- 1:300
autre <- int[-r3]
autre <- autre[-r2]
autre <- autre[-r1]
r <- c(r1, r2, r3)
vrai <- NULL
vrai[autre] <- FALSE
vrai[r] <- TRUE
vrai_intervals <- rle(vrai)

### lasso result (init)
plot(grid, vrai, type = "n", ylim = c(-1,1), xlab = "", ylab = "", axes = FALSE,
     xlim = c(0, 1.2))
axis(1)
est_intervals <- rle(res_sparse$sEDR[[1]][ ,1] != 0)
x0 <- x1 <- y0 <- y1 <- NULL
x0b <- x1b <- y0b <- y1b <- NULL
for (num in seq_along(vrai_intervals$lengths)) {
  if (vrai_intervals$values[num]) {
    if (num == 1) {
      x0 <- 1
      x1 <- vrai_intervals$lengths[num]
      y0 <- y1 <- 1
    } else {
      x0 <- c(x0, cumsum(vrai_intervals$lengths)[num-1])
      x1 <- c(x1, cumsum(vrai_intervals$lengths)[num])
      y0 <- c(y0, 1)
      y1 <- c(y1, 1)
    }
  }
}
for (num in seq_along(est_intervals$lengths)) {
  if (est_intervals$values[num]) {
    if (num == 1) {
      x0b <- -1
      x1b <- est_intervals$lengths[num]
      y0b <- y1b <- -1
    } else {
      x0b <- c(x0b, cumsum(est_intervals$lengths)[num-1])
      x1b <- c(x1b, cumsum(est_intervals$lengths)[num])
      y0b <- c(y0b, -1)
      y1b <- c(y1b, -1)
    }
  }
}
for (num in seq_along(x0)) {
  rect(grid[x0[num]], 0, grid[x1[num]], 1, col = "grey", density = NULL, 
       border = NA)
}
for (num in seq_along(x0b)) {
  rect(grid[x0b[num]], -1, grid[x1b[num]], 0, col = "pink", density = NULL, 
       border = NA)
}
segments(grid[x0], y0, grid[x1], y1, lwd=2)
segments(grid[x0b], y0b, grid[x1b], y1b, lwd = 2, col = "red")
lines(grid, rep(0, length(grid)))
# dev.print(png, file="../results/m2_intervals_init.png", width=600, height=400, res=150)

### best model
plot(grid, vrai, type="n", ylim=c(-1,1), xlab = "", ylab = "", axes = FALSE,
     xlim = c(0, 1.2))
axis(1)
est_intervals <- rle(res_sparse$sEDR[[bestCV]][ ,1] != 0)
x0 <- x1 <- y0 <- y1 <- NULL
x0b <- x1b <- y0b <- y1b <- NULL
for (num in seq_along(vrai_intervals$lengths)) {
  if (vrai_intervals$values[num]) {
    if (num == 1) {
      x0 <- 1
      x1 <- vrai_intervals$lengths[num]
      y0 <- y1 <- 1
    } else {
      x0 <- c(x0, cumsum(vrai_intervals$lengths)[num-1])
      x1 <- c(x1, cumsum(vrai_intervals$lengths)[num])
      y0 <- c(y0, 1)
      y1 <- c(y1, 1)
    }
  }
}
for (num in seq_along(est_intervals$lengths)) {
  if (est_intervals$values[num]) {
    if (num == 1) {
      x0b <- 1
      x1b <- est_intervals$lengths[num]
      y0b <- y1b <- -1
    } else {
      x0b <- c(x0b, cumsum(est_intervals$lengths)[num-1])
      x1b <- c(x1b, cumsum(est_intervals$lengths)[num])
      y0b <- c(y0b, -1)
      y1b <- c(y1b, -1)
    }
  }
}
for (num in seq_along(x0)) {
  rect(grid[x0[num]], 0, grid[x1[num]], 1, col = "grey", density = NULL, 
       border = NA)
}
for (num in seq_along(x0b)) {
  rect(grid[x0b[num]], -1, grid[x1b[num]], 0, col = "pink", density = NULL, 
       border = NA)
}
segments(grid[x0], y0, grid[x1], y1, lwd = 2)
segments(grid[x0b], y0b, grid[x1b], y1b, lwd = 2, col = "red")
lines(grid, rep(0, length(grid)))
# dev.print(png, file = "../results/m2_intervals.png", width = 600, height = 400, res = 150)