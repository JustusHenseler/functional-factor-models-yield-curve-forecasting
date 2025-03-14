source('R/fcst.RW.R')

oos.results.RW <- function(yield_data, max.h = 12, horizons = c(1,3,6,12), maturities = c(1,3,6,12,24,36,60,120,240,360)){
  
  full.densedata <- yield_data[[length(yield_data)]]$densedata
  
  start.training <- 120
  end.training <- dim(full.densedata)[1]-1  # T-h, i.e. the last object to base estimation on dependent on h
  period = start.training:end.training # periods to forecast consecutively
  
  ### Time Indices ####
  time.index <- seq(from = time(full.densedata)[1], by = 1/12, length.out = end.training + max.h) # extend time index to allow for forecasts out of sample window
  
  start.index <- time(full.densedata)[1]
  start.training.index <- time.index[start.training]
  end.training.index <- time.index[length(time.index)-1]
  
  ### Grid ####
  workinggrid <- yield_data[[length(yield_data)]]$workinggrid
  observationgrid <- yield_data[[length(yield_data)]]$observationgrid
  maturities_match <- match(maturities, workinggrid)
  
  len_period = length(period)
  len_grid = length(workinggrid)
  
  curve.extended <- ts(rbind(full.densedata, matrix(NA, nrow = max.h, ncol = length(workinggrid))), start=start.index, end = rev(time.index)[1], frequency = 12)
  
  empty.df <- data.frame(matrix(NA, nrow = len_period, ncol = len_grid))
  colnames(empty.df) <- workinggrid
  predicted.curves <- replicate(max.h, empty.df, simplify = FALSE)
  
  insample.fit <- data.frame(matrix(NA, nrow = len_period, ncol = 3))
  colnames(insample.fit) <- c("MSE", "RMSE", "MSE.L2")
  
  for (i in period){
    
    index <- which(period == i)
    cat('\rPeriod:', index, '/', len_period)
    
    pred.result = fcst.RW(max.h = max.h, obj = yield_data[[456]], end = i)
    
    insample.fit[index,"MSE"] <- pred.result$curve.mse
    insample.fit[index,"RMSE"] <- pred.result$curve.rmse
    insample.fit[index,"MSE.L2"] <- pred.result$curve.mse.l2
    
    for (h in 1:max.h){
      predicted.curves[[h]][index,] <- pred.result$predicted.curve[h,]
    }
  }
  
  cat(" - finished\n")
  
  insample.fit.avg <- colMeans(insample.fit, na.rm = TRUE)
  insample.fit.avg[1] <- sqrt(insample.fit.avg[1])
  names(insample.fit.avg) <- c("Root Avg. MSE", "Avg. RMSE", "Avg. MSE.L2")
  
  # turn predictions data frame into ts object
  predicted.curves <- lapply(1:max.h, function(h){
    predicted.curves <- ts(predicted.curves[[h]], start = time.index[start.training + h], end = time.index[end.training + h], frequency = 12)
    colnames(predicted.curves) <- workinggrid
    return(predicted.curves)
  })
  
  # calculate errors
  predicted.error.workinggrid <- lapply(1:max.h, function(h){
    error <- predicted.curves[[h]] - window(curve.extended, start = time.index[start.training + h], end = time.index[end.training + h])
    colnames(error) <- workinggrid
    return(error)
  })
  predicted.error.observationgrid <- lapply(1:max.h, function(h){
    return(predicted.error.workinggrid[[h]][,match(observationgrid, workinggrid)])
  })
  RMSFE.maturities <- lapply(1:max.h, function(h) sqrt(colMeans(predicted.error.observationgrid[[h]]^2, na.rm=TRUE)))
  RMSFE.overall <- unlist(lapply(1:max.h, function(h){
    return(sqrt(mean(colMeans(predicted.error.observationgrid[[h]]^2, na.rm=TRUE))))
  }))
  names(RMSFE.overall) <- paste0("h = ",1:h)
  
  output.overall <- data.frame(
    RMSFE = RMSFE.overall[horizons]
  )
  
  output.RMSFE.maturities <- do.call(rbind, lapply(RMSFE.maturities[horizons], `[`, maturities))
  
  rownames(output.RMSFE.maturities) <-  paste0("h = ", horizons)
  
  args <- as.character(pred.result$args)
  # specification <- paste0("K = ", args[1], ", p = ", args[2], ", AR = ", args[3], ", Macro = ",args[4], ", GARCH = ",args[5])
  
  show.output <- function() {    
    cat("Parameters: ", "RW", "\n")
    cat("\n### Overall RMSFE ###\n")
    print(output.overall)
    cat("\n\n### RMSFE over Maturities ###\n")
    print(output.RMSFE.maturities)
  }
  
  # show.output()
  
  output = list("args" = args,
                "predicted.error.workinggrid" = predicted.error.workinggrid,
                "predicted.error.observationgrid" = predicted.error.observationgrid,
                "predicted.curves" = predicted.curves,
                "insample.fit" = insample.fit,
                "insample.fit.avg" = insample.fit.avg,
                "RMSFE.overall" = RMSFE.overall,
                "RMSFE.maturities" = RMSFE.maturities,
                "output.overall" = output.overall,
                "output.RMSFE.maturities" = output.RMSFE.maturities,
                "show.output" = show.output()
  )
  
  return(output)
}
