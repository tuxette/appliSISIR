library(SISIR)
################################################################################
## import data
X <- as.matrix(read.table(file = "../data/climates.txt"))
Y <- as.matrix(read.table(file = "../data/yields.txt"))

# Actual relevant data: week 16 to week 41 => days112 to 287
# Everything before p1 or after p2 should be zero
p1 <- 0.27
p2 <- 0.84

p <- dim(X)[2]
n <- length(Y)
tsteps <- seq(0, 1, length.out = p)

# slices
H <- 10
slices_bounds <- quantile(Y, probs = seq(0, 1, length = H+1))
y_class <- cut(Y, breaks = slices_bounds, labels = FALSE, include.lowest = TRUE)
colH <- rainbow(length(slices_bounds))[y_class]
keep <- which(!duplicated(y_class))

# plots
par(mfrow = c(1,1), mar=c(4,4,1,1))
plot(NA, type = "l", xlim = range(tsteps), ylim = range(X), xlab = "time",
     ylab = "Evapotranspiration")
for (ind in keep) {
  lines(tsteps, X[ind, ], col = colH[ind])
}
abline(v = c(p1, p2), lty = 4)
# dev.off(png, file="../results/sunflo-data.png", width=800, height=500, res=150)

################################################################################
## Ridge step  
set.seed(150)
res_ridge <- ridgeSIR(X, Y, H = 10, d = 2, mu2 = 1000)

################################################################################
## Sparse estimation
res_sparse <- SISIR(res_ridge, parallel = FALSE)

################################################################################
## Analysis of results
bestCV <- which.min(res_sparse$quality$CV + 1e12*(res_sparse$quality$nbint>200))

################################################################################
## Plots
# Criteria evolution
par(mfrow=c(1,1))
plot(res_sparse$quality$nbint, res_sparse$quality$CVerror, col="red", type="l",
     xlab="Number of intervals", ylab="CV error")
abline(v=res_sparse$quality$nbint[bestCV])
# dev.print(png, file="../results/sunflo-cv.png", width=1000, height=500, res=150)

# lasso and best models
par(mar = c(4,2,2,2))
grid <- seq(0, 1, length = p)

### lasso result (init)
est_intervals <- rle(res_sparse$sEDR[[1]][ ,1] != 0)
x0b <- x1b <- y0b <- y1b <- NULL
for (num in seq_along(est_intervals$lengths)) {
  if (est_intervals$values[num]) {
    if (num == 1) {
      x0b <- -1
      x1b <- est_intervals$lengths[num]
      y0b <- y1b <- 1
    } else {
      x0b <- c(x0b, cumsum(est_intervals$lengths)[num-1])
      x1b <- c(x1b, cumsum(est_intervals$lengths)[num])
      y0b <- c(y0b, 1)
      y1b <- c(y1b, 1)
    }
  }
}
plot(NA, type = "n", ylim = c(-1,1), xlab = "", ylab = "", axes = FALSE,
     xlim = c(0, 1))
axis(1)
for (num in seq_along(x0b)) {
  rect(grid[x0b[num]], 1, grid[x1b[num]], 0, col = "lightblue", density = NULL, 
       border = NA)
}
segments(grid[x0b], y0b, grid[x1b], y1b, lwd = 2, col = "darkblue")
lines(grid, rep(0, length(grid)))

### best model
est_intervals <- rle(res_sparse$sEDR[[bestCV]][ ,1] != 0)
x0b <- x1b <- y0b <- y1b <- NULL
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
for (num in seq_along(x0b)) {
  rect(grid[x0b[num]], -1, grid[x1b[num]], 0, col = "pink", density = NULL, 
       border = NA)
}
segments(grid[x0b], y0b, grid[x1b], y1b, lwd = 2, col = "red")
lines(grid, rep(0, length(grid)))

abline(v = c(p1,p2), lty = 4, lwd = 2)
# dev.print(png, file="../results/sunflo-models.png", width=1000, height=1000, res=150)
