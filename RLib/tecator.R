library(SISIR)

## PACKAGE
library(e1071)

##### Load data ---------------------------------------------------------------
tecator_init <- as.matrix(read.table("../data/tecator.txt"))

tecator <- matrix(0, ncol = 125, nrow = 240)
for (indiv in 1:240) {
  for (size in 1:25) {
    indcol <- (((size-1)*5+1):(size*5))
    tecator[indiv,indcol] <- tecator_init[(indiv-1)*25+size, ]
  }
}
spectra <- tecator[ ,1:100]
fat <- tecator[ ,124]
tsteps <- seq(850, 1050, length = 100)

# derivatives
sp.d1 <- (spectra[ ,2:100] - spectra[ ,1:99])/(tsteps[2]-tsteps[1])
sp.d2 <- (sp.d1[ ,2:99] - sp.d1[ ,1:98])/(tsteps[2]-tsteps[1])

# prepare data
set.seed(17031058)
n <- nrow(spectra)
ntrain <- 200
ntest <- n - ntrain
sel_dat <- "1"
x <- switch(sel_dat,
            "0" = spectra,
            "1" = sp.d1,
            "2" = sp.d2)
y <- fat
active_steps <- tsteps[1:ncol(x)]

# slices
H <- 10
slices_bounds <- quantile(y, probs = seq(0, 1, length=H+1))
y_class <- cut(y, breaks=slices_bounds, labels = FALSE, include.lowest = TRUE)
colH <- rainbow(length(slices_bounds))[y_class]

# plots
op <- par()
par(mfrow = c(1,2))
plot(tsteps, spectra[1, ], type="l", ylim = range(spectra), xlab = "wavelength",
     ylab = "absorbance", col = colH[1], main = "original spectra")
for (ind in 2:nrow(spectra)) {
  lines(tsteps, spectra[ind, ], col = colH[ind])
}
plot(tsteps[1:ncol(sp.d1)], sp.d1[1, ], type="l", ylim = range(sp.d1), 
     xlab = "wavelength", ylab = "absorbance", col = colH[1], 
     main = "derivatives")
for (ind in 2:nrow(spectra)) {
  lines(tsteps[1:ncol(sp.d1)], sp.d1[ind, ], col = colH[ind])
}
# dev.print(png, file = "../results/tecator-data.png", width = 1000, height = 500)

##### Process data once -------------------------------------------------------
# select mu2 and d
list_mu2 <- 10^(-6:0)
listH <- c(H)
list_d <- 1:4
set.seed(1603)
res_tune <- tune.ridgeSIR(x, y, listH, list_mu2, list_d, nfolds = 10,
                          parallel = FALSE)

printable_res <- matrix(res_tune[ ,"cverror"], ncol = length(list_d),
                        byrow = TRUE)
rownames(printable_res) <- log10(list_mu2)
colnames(printable_res) <- list_d
par(op)
plot(log10(list_mu2), printable_res[ ,ncol(printable_res)], main = "",
     ylab = "CV error", xlab = expression(log[10](mu[2])), type = "b",
     pch = "+")
chosen_mu <- which.min(printable_res[ ,ncol(printable_res)])
## chosen_mu: 10(-4)

## selection of d
printable_res2 <- matrix(res_tune[res_tune$H == 10, "Rd"], 
                         ncol = length(list_d), byrow = TRUE)
rownames(printable_res2) <- log10(list_mu2)
colnames(printable_res2) <- list_d
plot(list_d, printable_res2[chosen_mu, ], type = "b", pch = "+", xlab = "d",
     ylab = expression(R(d)))
d <- 1

## final check
plot(log10(list_mu2), printable_res[ ,d], main = "", ylab = "CV error",
     xlab = expression(log[10](mu[2])), type = "b", pch = "+")
which.min(printable_res[ ,ncol(printable_res)]) ## consistent


##### Load data ---------------------------------------------------------------
set.seed(22111723)
seeds <- sample(1:(10^5), 100, replace = FALSE)
registerDoParallel(cores = detectCores() - 1)
res <- foreach(repet = 1:100, .combine = c) %dopar% {
  set.seed(seeds[repet])
  ind_test <- sample(1:n, ntest, replace = FALSE)
  # train/test split
  x_train <- x[- ind_test, ]
  y_train <- y[- ind_test]
  x_test <- x[ind_test,]
  y_test <- y[ind_test]

  ##### Ridge -----------------------------------------------------------------
  print("ridge...")
  res_ridge <- ridgeSIR(x_train, y_train, H=H, d=d, mu2=list_mu2[chosen_mu])

  ##### Lasso -------------------------------------------------------------------
  ## Search for relevant intervals
  print("sparse...")
  res_sparse <- SISIR(res_ridge, parallel = FALSE, minint = 20)

  ind_opt <- which.min(res_sparse$quality$CVerror)
  print("regression...")
  beta1 <- res_sparse$sEDR[[ind_opt]][ ,1]
  preds <- x_test %*% beta1
  res_tune <- tune.svm(preds, y_test, gamma = 10^(-3:3), cost = 10^(-3:3))
  
  pred <- res_tune$best.model$fitted
  pred <- mean((pred - y_test)^2)
  
  return(pred)
}

summary(res)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.090   4.087   4.913   5.536   6.206  14.050

# save(list=ls(), file = "../results/result_tecator.rda")

##### Analyze one result ------------------------------------------------------
set.seed(seeds[3])
ind_test <- sample(1:n, ntest, replace = FALSE)
# train/test split
x_train <- x[- ind_test, ]
y_train <- y[- ind_test]
x_test <- x[ind_test,]
y_test <- y[ind_test]

# ridge
res_ridge <- ridgeSIR(x_train, y_train, H=H, d=d, mu2=list_mu2[chosen_mu])

# sparse
res_sparse <- SISIR(res_ridge, parallel = FALSE, minint = 20)

ind_opt <- which.min(res_sparse$quality$CVerror)
# print(res_sparse$quality$nbint[ind_best]): 89
hatA <- as.vector(res_sparse$sEDR[[ind_best]])
nt <- length(hatA)
success <- hatA != 0
starts <- !success[2:(nt-1)] & success[3:nt]
starts <- c(success[1], starts, FALSE)
ends <- !success[2:nt] & success[1:(nt-1)]
ends <- c(FALSE, ends)
plot(tsteps[1:ncol(sp.d1)], sp.d1[1, ], type="n", ylim = range(sp.d1), 
     xlab = "wavelength", ylab = "absorbance", main = "derivatives")
rect(active_steps[starts], rep(min(sp.d1), length(sp.d1)), 
     active_steps[ends], rep(max(sp.d1), length(sp.d1)), col = "grey",
     density = NULL, border = NA)
for (ind in 1:nrow(spectra)) {
  lines(tsteps[1:ncol(sp.d1)], sp.d1[ind, ], col = colH[ind])
}
# dev.print(png, file = "../results/tecator-relevantint.png", width = 500, height = 500)

hatA <- as.vector(res_sparse$sEDR[[1]])
nt <- length(hatA)
success <- hatA != 0
starts <- !success[2:(nt-1)] & success[3:nt]
starts <- c(success[1], starts, FALSE)
ends <- !success[2:nt] & success[1:(nt-1)]
ends <- c(FALSE, ends)
plot(tsteps[1:ncol(sp.d1)], sp.d1[1, ], type="n", ylim = range(sp.d1), 
     xlab = "wavelength", ylab = "absorbance", main = "derivatives")
rect(active_steps[starts], rep(min(sp.d1), length(sp.d1)), 
     active_steps[ends], rep(max(sp.d1), length(sp.d1)), col = "grey",
     density = NULL, border = NA)
for (ind in 1:nrow(spectra)) {
  lines(tsteps[1:ncol(sp.d1)], sp.d1[ind, ], col = colH[ind])
}
# dev.print(png, file = "../results/tecator-relevantint-init.png", width = 500, height = 500)