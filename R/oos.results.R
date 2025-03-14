library("vars")

oos.results <- function(model, yield_data, macro_factors, max.h = 12, horizons = c(1,3,6,12), maturities = c(1,3,6,12,24,36,60,120,240,360), IC, IC_type = "BIC", GARCH = TRUE){

  if (!missing(IC) && !IC_type %in% c("BIC", "HQC")) {
    stop("IC_type can only be equal to 'BIC' oder 'HQC'.")
  }
  
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
  predicted.curves <- predicted.curve.FIs.lower <- predicted.curve.FIs.upper <- replicate(max.h, empty.df, simplify = FALSE)
  
  if (GARCH){
    predicted.curve.FIs.GARCH.lower <- predicted.curve.FIs.GARCH.upper <- replicate(max.h, empty.df, simplify = FALSE)
  }
  
  insample.fit <- data.frame(matrix(NA, nrow = len_period, ncol = 3))
  colnames(insample.fit) <- c("MSE", "RMSE", "MSE.L2")
  
  for (i in period){
    
    index <- which(period == i)
    cat('\rPeriod:', index, '/', len_period)
    
    args <- c(max.h = max.h, obj = yield_data[i], model$args, end = i)
    
    if (!missing(macro_factors)) {
      args$macro_data <- macro_factors[[i]]  
    }
    if (!missing(IC) && IC_type == "BIC") {
      args$K <- IC[[i]][1,1]    
      args$p <- IC[[i]][1,2]    
    }
    if (!missing(IC) && IC_type == "HQC") {
      args$K <- IC[[i]][2,1]    
      args$p <- IC[[i]][2,2]    
    }
    if (GARCH) {
      args$GARCH <- TRUE   
    }
    
    pred.result = do.call(model$func, args)
    
    insample.fit[index,"MSE"] <- pred.result$curve.mse
    insample.fit[index,"RMSE"] <- pred.result$curve.rmse
    insample.fit[index,"MSE.L2"] <- pred.result$curve.mse.l2
    
    for (h in 1:max.h){
      predicted.curves[[h]][index,] <- pred.result$predicted.curve[h,]
      predicted.curve.FIs.lower[[h]][index,] <- pred.result$predicted.curve.FI.lower[h,]
      predicted.curve.FIs.upper[[h]][index,] <- pred.result$predicted.curve.FI.upper[h,]
      if (GARCH){
      predicted.curve.FIs.GARCH.lower[[h]][index,] <- pred.result$predicted.curve.FI.GARCH.lower[h,]
      predicted.curve.FIs.GARCH.upper[[h]][index,] <- pred.result$predicted.curve.FI.GARCH.upper[h,]
      }
    }
  }
  
  cat(" - finished\n")
  args <- as.character(pred.result$args)
  
  insample.fit.avg <- colMeans(insample.fit, na.rm = TRUE)
  insample.fit.avg[1] <- sqrt(insample.fit.avg[1])
  names(insample.fit.avg) <- c("Root Avg. MSE", "Avg. RMSE", "Avg. MSE.L2")
  
  # turn predictions data frame into ts object
  predicted.curves <- lapply(1:max.h, function(h){
    predicted.curves <- ts(predicted.curves[[h]], start = time.index[start.training + h], end = time.index[end.training + h], frequency = 12)
    colnames(predicted.curves) <- workinggrid
    return(predicted.curves)
  })
  predicted.curve.FIs.lower <- lapply(1:max.h, function(h){
    predicted.curve.FIs.lower <- ts(predicted.curve.FIs.lower[[h]], start = time.index[start.training + h], end = time.index[end.training + h], frequency = 12)
    colnames(predicted.curve.FIs.lower) <- workinggrid
    return(predicted.curve.FIs.lower)
  })
  predicted.curve.FIs.upper <- lapply(1:max.h, function(h){
    predicted.curve.FIs.upper <- ts(predicted.curve.FIs.upper[[h]], start = time.index[start.training + h], end = time.index[end.training + h], frequency = 12)
    colnames(predicted.curve.FIs.upper) <- workinggrid
    return(predicted.curve.FIs.upper)
  })
  
  if (GARCH){
    predicted.curve.FIs.GARCH.lower <- lapply(1:max.h, function(h){
      predicted.curve.FIs.GARCH.lower <- ts(predicted.curve.FIs.GARCH.lower[[h]], start = time.index[start.training + h], end = time.index[end.training + h], frequency = 12)
      colnames(predicted.curve.FIs.GARCH.lower) <- workinggrid
      return(predicted.curve.FIs.GARCH.lower)
    })
    predicted.curve.FIs.GARCH.upper <- lapply(1:max.h, function(h){
      predicted.curve.FIs.GARCH.upper <- ts(predicted.curve.FIs.GARCH.upper[[h]], start = time.index[start.training + h], end = time.index[end.training + h], frequency = 12)
      colnames(predicted.curve.FIs.GARCH.upper) <- workinggrid
      return(predicted.curve.FIs.GARCH.upper)
    })
  }
  
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
  
  coverage <- lapply(1:max.h, function(h){
    curve.oos <- window(curve.extended, start = time.index[start.training + h], end = time.index[end.training + h])
    coverage <- ((curve.oos >= predicted.curve.FIs.lower[[h]]) & (curve.oos <= predicted.curve.FIs.upper[[h]]))
    colnames(coverage) <- workinggrid
    return(coverage)
  })
  
  if(GARCH){
    coverage.GARCH <- lapply(1:max.h, function(h){
      curve.oos <- window(curve.extended, start = time.index[start.training + h], end = time.index[end.training + h])
      coverage.GARCH <- ((curve.oos >= predicted.curve.FIs.GARCH.lower[[h]]) & (curve.oos <= predicted.curve.FIs.GARCH.upper[[h]]))
      colnames(coverage.GARCH) <- workinggrid
      return(coverage.GARCH)
    })
  }
  
  coverage.maturities <- lapply(coverage, function(x) colMeans(x, na.rm = TRUE))
  coverage.overall <-  unlist(lapply(coverage.maturities, function(x) mean(x, na.rm = TRUE))) ##
  names(coverage.overall) <- paste0("h = ", 1:max.h)
  
  if(GARCH){
    coverage.GARCH.maturities <- lapply(coverage.GARCH, function(x) colMeans(x, na.rm = TRUE))
    coverage.GARCH.overall <-  unlist(lapply(coverage.GARCH.maturities, function(x) mean(x, na.rm = TRUE))) ##
    names(coverage.GARCH.overall) <- paste0("h = ", 1:max.h)
  }
  
  output.overall <- data.frame(
    RMSFE = RMSFE.overall[horizons],
    Coverage = coverage.overall[horizons]
  )
  if (GARCH){
    output.overall["Coverage DCC-GARCH"] <- coverage.GARCH.overall[horizons]
  }
  
  output.RMSFE.maturities <- do.call(rbind, lapply(RMSFE.maturities[horizons], `[`, maturities))
  output.coverage.maturities <- do.call(rbind, lapply(coverage.maturities[horizons], `[`, maturities_match))
  rownames(output.RMSFE.maturities) <-  rownames(output.coverage.maturities) <- paste0("h = ", horizons)
  if(GARCH){
  output.coverage.GARCH.maturities <- do.call(rbind, lapply(coverage.GARCH.maturities[horizons], `[`, maturities_match))
  rownames(output.coverage.GARCH.maturities) <- paste0("h = ", horizons)
  }
  
  # colnames(output.coverage.maturities) <- paste0(maturities)
  
  specification <- paste0("K = ", args[1], ", p = ", args[2], ", AR = ", args[3], ", Macro = ",args[4], ", GARCH = ",args[5])
  
  show.output.fun <- function() {    
    cat("Parameters: ", specification, "\n")
    cat("\n### Overall RMSFE & Coverage ###\n")
    print(output.overall)
    cat("\n\n### RMSFE over Maturities ###\n")
    print(output.RMSFE.maturities)
    cat("\n\n### Coverage over Maturities ###\n")
    print(output.coverage.maturities)
  }
  
  if(GARCH){
  show.output.fun <- function() {    
    cat("Parameters: ", specification, "\n")
    cat("\n### Overall RMSFE & Coverage ###\n")
    print(output.overall)
    cat("\n\n### RMSFE over Maturities ###\n")
    print(output.RMSFE.maturities)
    cat("\n\n### Coverage over Maturities ###\n")
    print(output.coverage.maturities)
    cat("\n\n### Coverage over Maturities (DCC-GARCH) ###\n")
    print(output.coverage.GARCH.maturities)
  }
  }
  
  show.output <- capture.output(show.output.fun())
  
  # show.output()
  
  output = list("args" = args,
                "predicted.error.workinggrid" = predicted.error.workinggrid,
                "predicted.error.observationgrid" = predicted.error.observationgrid,
                "predicted.curves" = predicted.curves,
                "predicted.curve.FIs.lower" = predicted.curve.FIs.lower,
                "predicted.curve.FIs.upper" = predicted.curve.FIs.upper,
                "coverage" = coverage,
                "coverage.overall" = coverage.overall,
                "coverage.maturities" = coverage.maturities,
                "insample.fit" = insample.fit,
                "insample.fit.avg" = insample.fit.avg,
                "RMSFE.overall" = RMSFE.overall,
                "RMSFE.maturities" = RMSFE.maturities,
                "output.overall" = output.overall,
                "output.RMSFE.maturities" = output.RMSFE.maturities,
                "output.coverage.maturities" = output.coverage.maturities,
                "show.output" = show.output
  )
  
  if (GARCH){
    output$coverage.GARCH <- coverage.GARCH
    output$coverage.GARCH.overall <- coverage.GARCH.overall
    output$coverage.GARCH.maturities <- coverage.GARCH.maturities
    output$output.coverage.GARCH.maturities <- output.coverage.GARCH.maturities
    output$predicted.curve.FIs.GARCH.upper <- predicted.curve.FIs.GARCH.upper
    output$predicted.curve.FIs.GARCH.lower <- predicted.curve.FIs.GARCH.lower
  }
  
  
  return(output)
}
