library(xts)
library(vars)

load.FRED_MD <- function(){
  
  FRED.MD_data <- read.csv("data/FRED-MD_current.csv", skip = 0, row.names = 1, header =TRUE)
  
  transformation_factors <- FRED.MD_data[1,]
  FRED.MD_data <- FRED.MD_data[-1,]
  indices <- as.Date(rownames(FRED.MD_data), format = "%m/%d/%Y")
  rownames(FRED.MD_data) <- indices

  # Apply transformations as given by Appendix of McCracken and Ng (2015)
  FRED.MD_data <- as.data.frame(lapply(1:ncol(FRED.MD_data), function(i) transform_column(FRED.MD_data[[i]], transformation_factors[i])))
  
  # Remove strong outliers (z-score bigger than 4) and replace by NA
  FRED.MD_data <- as.data.frame(lapply(FRED.MD_data, remove_outliers_z, threshold = 4))
  
  # columns_na_percentage <- apply(is.na(FRED.MD_data), 2, mean) # percentage of NAs per column
  FRED.MD_data <- apply(FRED.MD_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
  # FRED.MD_data <- scale(FRED.MD_data)
  FRED.MD_data <- ts(FRED.MD_data,start=c(1986,1),end=c(2023,12),frequency=12)
  
  return(FRED.MD_data)
}

transform_column <- function(column, factor) {
  n <- length(column)
  
  if      (factor == 1) {return(column)}
  else if (factor == 2) {return(c(NA, diff(column)))} 
  else if (factor == 3) {return(c(NA, NA, diff(diff(column))))} 
  else if (factor == 4) {return(log(column))} 
  else if (factor == 5) {return(c(NA, diff(log(column))))} 
  else if (factor == 6) {return(c(NA, NA, diff(diff(log(column)))))} 
  else if (factor == 7) {return(c(NA, (column[-1] / column[-n]) - 1))} 
  else                  {stop("Factor invalid!")
  }
}

remove_outliers_z <- function(series, threshold = 3) {
  mean_val <- mean(series, na.rm = TRUE)
  sd_val <- sd(series, na.rm = TRUE)
  z_scores <- (series - mean_val) / sd_val
  
  series[abs(z_scores) > threshold] <- NA
  return(series)
}
