library(dffm)
library(vars)
library(rugarch)
library(rmgarch)

fcst.RW <- function(obj, 
                   start = 1, 
                   end = 120, 
                   max.h = 12){

densedata <- obj$densedata
workinggrid <- obj$workinggrid
observationgrid <- obj$observationgrid

### Time Indices ####
time.index <- time(densedata) # get time index of factors
time.index.extended <- seq(from = time.index[1], by = 1/12, length.out = length(time.index) + max.h) # extend time index to allow for forecasts out of sample window

start.training <- time.index[start]
end.training <- time.index[end]
start.fitted <- time.index[start+1]

start.forecast <- time.index.extended[end+1]
end.forecast <- time.index.extended[end+max.h]

# RW
curve.fitted <- lag(densedata,1)[1:dim(densedata)[1]-1,]
curve.fitted <- ts(curve.fitted, start = start.fitted, end=end.training, frequency = 12)

curve.residuals <- curve.fitted - window(densedata, start = start.fitted, end=end.training, frequency = 12)
# curve.residuals <- diff(densedata,1)
curve.residuals.observationgrid <- curve.residuals[,match(observationgrid, workinggrid)]

binwidth <- (rev(workinggrid)[1]-(workinggrid[1]))/length(workinggrid)
curve.mse.l2 <- mean(rowSums(curve.residuals^2)*binwidth)
curve.mse <- mean(curve.residuals.observationgrid^2, na.rm=TRUE)
curve.rmse <- sqrt(curve.mse)

last_curve <- window(densedata, start = end.training, end = end.training)
predicted.curve <- do.call(rbind, replicate(max.h, last_curve , simplify = FALSE))
predicted.curve <- ts(predicted.curve, start = start.forecast, end=end.forecast, frequency = 12)

args <- list('start' = start, 'end' = end, 'max.h' = max.h)

output = list("args" = args,
              "curve.fitted" = curve.fitted,
              "curve.residuals" = curve.residuals,
              "curve.mse.l2" = curve.mse.l2,
              "curve.mse" = curve.mse,
              "curve.rmse" = curve.rmse,
              "predicted.curve" = predicted.curve
)

return(output)
  
}