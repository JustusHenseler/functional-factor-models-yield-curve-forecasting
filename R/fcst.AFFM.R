library(dffm)
library(vars)
library(rugarch)
library(rmgarch)

fcst.AFFM = function(obj, 
                        macro_data,
                        GARCH = FALSE,
                        K = 3, 
                        p = 1, 
                        AR = FALSE, 
                        start = 1, 
                        end = 120, 
                        max.h = 12, 
                        restrict_macro = FALSE){
  
macro <- !missing(macro_data) # check if macro data is to be incorporated in forecast

densedata <- obj$densedata
loadings <- obj$eigenfunctions[,1:K,drop=F]
factors <- obj$scores.centered[,1:K,drop=F] # determine loadings & factors

workinggrid <- obj$workinggrid # select workinggrid from fdaobj
observationgrid <- obj$observationgrid # select observationgrid from fdaobj

### Time Indices ####
time.index <- time(factors) # get time index of factors
time.index.extended <- seq(from = time.index[1], by = 1/12, length.out = length(time.index) + max.h) 
    # extend time index to allow for forecasts out of sample window 

start.training <- time.index[start]
end.training <- time.index[end]
start.fitted <- time.index[start+p]

start.forecast <- time.index.extended[end+1]
end.forecast <- time.index.extended[end+max.h]

### Select factors for training data ####
densedata <- window(densedata, start = start.training, end = end.training)
factors <- window(factors, start = start.training, end = end.training)
colnames.factors <- paste0("F",1:K)
colnames(factors) <- colnames.factors

### Calculate curve residuals
curve.epsilon <- t(replicate(dim(factors)[1], obj$meanfunction)) + t(loadings %*% t(factors)) - densedata
curve.epsilon.sample.variances <- apply(curve.epsilon, 2, var)

### Add macro factors to training data if macro data is given ####
if (macro){
  no.macro <- ncol(macro_data)
  colnames.macro_data <- paste0("M",1:no.macro)
  colnames(macro_data) <- colnames.macro_data
  macro_data <- window(macro_data, start = start.training, end = end.training)
  data <- cbind(factors, macro_data)
  colnames.data <-  c(colnames.factors, colnames.macro_data)
  colnames(data) <- colnames.data
  restriction_AR_macro_one_lag <- if (restrict_macro) rbind(cbind(diag(K), matrix(1, K, no.macro)), cbind(matrix(0, no.macro, K), diag(no.macro))) else rbind(cbind(diag(K), matrix(1, K, no.macro)), cbind(matrix(1, no.macro, K),matrix(1, no.macro, no.macro)))
  restriction_matrix_AR_macro <- do.call(cbind, replicate(p, restriction_AR_macro_one_lag, simplify = FALSE))
} else {
  data <- factors
  colnames.data <-  colnames.factors
  colnames(data) <- colnames.data
  restriction_matrix_AR <- do.call(cbind, replicate(p, diag(ncol(data)), simplify = FALSE))
}

### Estimate model ####
var_model <- VAR(data, p = p, type = "none")
var_model <- if (AR & !macro) restrict(var_model, method = "manual", resmat = restriction_matrix_AR) else var_model
var_model <- if (AR & macro) restrict(var_model, method = "manual", resmat = restriction_matrix_AR_macro) else var_model

### In-sample metrics of data ####
data.fitted <- fitted(var_model)
data.residuals <- residuals(var_model)
data.COV <- summary(var_model)$covres
  
### In-sample metrics of approximate factors ####
factors.fitted <- data.fitted[,colnames.factors]
factors.residuals <- data.residuals[,colnames.factors] 
# factors.COV <- ((nrow(factors.residuals) - 1) / (nrow(factors.residuals) - p * K)) * cov(factors.residuals)
factors.COV <- data.COV[colnames.factors,colnames.factors]

# Calculate curve residuals on workinggrid & observationgrid to calculate evaluation metrics
length.fitted <- dim(factors.fitted)[1]
curve.insample <- window(densedata, start = start.fitted, end=end.training)
curve.fitted <- t(replicate(length.fitted, obj$meanfunction)) + factors.fitted %*% t(loadings)
curve.fitted <- ts(curve.fitted, start = start.fitted, end=end.training, frequency = 12)
curve.residuals <- curve.insample - curve.fitted
curve.residuals.observationgrid <- curve.residuals[,match(observationgrid, workinggrid)]

binwidth = (rev(workinggrid)[1]-(workinggrid[1]))/length(workinggrid)
curve.mse.l2 = mean(rowSums(curve.residuals^2)*binwidth)
curve.mse = mean(curve.residuals.observationgrid^2, na.rm=TRUE)
curve.rmse = sqrt(curve.mse)

### Predict factor data
var_forecast <- predict(var_model, n.ahead = max.h)
predictions <- var_forecast$fcst
predicted.data <- data.frame(lapply(predictions, function(x) x[, "fcst"]))
predicted.factors <- predicted.data[colnames.factors]

### Predict covariance of forecast (given Wolds MA representation and forecasted error covariance matrices) ####
# predicted.error.COV <- array(replicate(max.h, data.COV), dim = c(nrow(data.COV), ncol(data.COV), max.h))
# predicted.error.COV.reversed <- predicted.error.COV[,,rev(1:max.h)] # Reverse order of error covariance for easier computations
data.Phi <- Phi(var_model, nstep = max.h) # Calculated Coefficients of Wold MA Representation
Wold.summands <- lapply(c(1:max.h), function(i) data.Phi[, ,i] %*% data.COV %*% t(data.Phi[, ,i])) # Calculate summands of Wold MA representattion
predicted.data.COV <- lapply(c(1:max.h), function(i) {Reduce("+", Wold.summands[1:i])}) # Save forecasted variance for each h in list
  # names(predicted.data.COV) <- paste0("h = ",1:max.h)
predicted.data.COV <- lapply(predicted.data.COV, function(x) {
  rownames(x) <- colnames(x) <- colnames.data
  return(x)
})
predicted.factors.COV <- lapply(predicted.data.COV, function(x) {
  x <- x[colnames.factors, colnames.factors]
  return(x)
}) # only consider forecasted covariance matrix of approximate factors

### Forecast covariance matrix via GARCH ####
if (GARCH){
  garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  dcc_spec <- dccspec(
    uspec = multispec(replicate(ncol(data), garch_spec)),
    dccOrder = c(1, 1),
    distribution = "mvnorm"
  )
  dcc_fit <- dccfit(dcc_spec, data = data.residuals)
  dcc_forecast <- dccforecast(dcc_fit, n.ahead = max.h)
  
  predicted.error.COV.GARCH <- dcc_forecast@mforecast[["H"]][[1]]
  
  calculate_predicted.variance <- function(h){
    summands <- list()
    for (i in 1:h) {
      summands[[i]] <- data.Phi[, ,rev(1:h)[i]] %*% predicted.error.COV.GARCH[, ,i] %*% t(data.Phi[, ,rev(1:h)[i]])
    }
    return(Reduce("+", summands))
  }
  predicted.data.COV.GARCH <- lapply(c(1:max.h), calculate_predicted.variance) # Save forecasted variance for each h in list
  # names(predicted.data.COV.GARCH) <- paste0("h = ",1:max.h)
  predicted.data.COV.GARCH <- lapply(predicted.data.COV.GARCH, function(x) {
    rownames(x) <- colnames(x) <- colnames.data
    return(x)
  })
  predicted.factors.COV.GARCH <- lapply(predicted.data.COV.GARCH, function(x) {
    x <- x[colnames.factors, colnames.factors]
    return(x)
  }) # only consider forecasted covariance matrix of approximate factors
} 



#### Predict point forecast and FI of curves for horizons 1:max.h
predicted.curve <- data.frame(matrix(NA, nrow = max.h, ncol = dim(densedata)[2]))
colnames(predicted.curve) <- workinggrid

predicted.curve.variance <- predicted.curve
predicted.curve.FI.lower <- predicted.curve
predicted.curve.FI.upper <- predicted.curve

if (GARCH){
  predicted.curve.variance.GARCH <- predicted.curve
  predicted.curve.FI.GARCH.lower <- predicted.curve
  predicted.curve.FI.GARCH.upper <- predicted.curve
}

for (i in 1:max.h){
  predicted.curve[i,] <- as.numeric(predicted.factors[i,]) %*% t(loadings) + obj$meanfunction
  # predicted.curve[i,] <- t(loadings %*% predicted.factors.current) + fdaobj$meanfunction
  
  predicted.curve.variance[i,] <- diag(loadings %*% predicted.factors.COV[[i]]  %*% t(loadings)) + curve.epsilon.sample.variances
  predicted.curve.FI.lower[i,] <- predicted.curve[i,] - 1.959964*sqrt(predicted.curve.variance[i,])
  predicted.curve.FI.upper[i,] <- predicted.curve[i,] + 1.959964*sqrt(predicted.curve.variance[i,])
  
  if (GARCH){
    predicted.curve.variance.GARCH[i,] <- diag(loadings %*% predicted.factors.COV.GARCH[[i]]  %*% t(loadings)) + curve.epsilon.sample.variances # CHECKKKKKKKK !!!!!!!!
    predicted.curve.FI.GARCH.lower[i,] <- predicted.curve[i,] - 1.959964*sqrt(predicted.curve.variance.GARCH[i,])
    predicted.curve.FI.GARCH.upper[i,] <- predicted.curve[i,] + 1.959964*sqrt(predicted.curve.variance.GARCH[i,])
  }
}


### Creat ts objects
predicted.factors <- ts(predicted.factors, start = start.forecast, end=end.forecast, frequency = 12)
predicted.curve <- ts(predicted.curve, start = start.forecast, end=end.forecast, frequency = 12)
predicted.curve.FI.lower <- ts(predicted.curve.FI.lower, start = start.forecast, end=end.forecast, frequency = 12)
predicted.curve.FI.upper <- ts(predicted.curve.FI.upper, start = start.forecast, end=end.forecast, frequency = 12)


args <- list('K' = K, 'p' = p, 'AR' = AR, 'macro' = macro, 'GARCH' = GARCH, 'start' = start, 'end' = end, 'max.h' = max.h)

output <- list("args" = args,
           "factors.fitted" = factors.fitted,
           "factors.residuals" = factors.residuals,
           "curve.fitted" = curve.fitted,
           "curve.residuals" = curve.residuals,
           "curve.mse.l2" = curve.mse.l2,
           "curve.mse" = curve.mse,
           "curve.rmse" = curve.rmse,
           "predicted.factors" = predicted.factors,
           "precicted.factors.COV" = predicted.factors.COV,
           "predicted.curve" = predicted.curve,
           "predicted.curve.FI.upper" = predicted.curve.FI.upper,
           "predicted.curve.FI.lower" = predicted.curve.FI.lower
)

if (GARCH){
  predicted.curve.FI.GARCH.lower <- ts(predicted.curve.FI.GARCH.lower, start = start.forecast, end=end.forecast, frequency = 12)
  predicted.curve.FI.GARCH.upper <- ts(predicted.curve.FI.GARCH.upper, start = start.forecast, end=end.forecast, frequency = 12)
  output$predicted.curve.FI.GARCH.upper <- predicted.curve.FI.GARCH.upper
  output$predicted.curve.FI.GARCH.lower <- predicted.curve.FI.GARCH.lower
}

return(output)
}