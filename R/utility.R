factor_loading_2 <- function(tau, lambda){return((1-exp(-lambda*tau))/(lambda*tau))}
factor_loading_3 <- function(tau, lambda){return(((1-exp(-lambda*tau))/(lambda*tau))-exp(-lambda*tau))}

DNS_factors <- function(data, lambda = 0.0609, maturities){
  
  factor_loading_2_vec <- sapply(maturities/12, function(i) factor_loading_2(i, lambda))
  factor_loading_3_vec <- sapply(maturities/12, function(i) factor_loading_3(i, lambda))
  
  return(t(apply(data, 1, function(i) lm(as.numeric(i) ~  factor_loading_2_vec + factor_loading_3_vec)$coefficients)))
}

descriptiveStats <- function(data){
  
  df <- data.frame(   colMeans(data),
                      apply(data, 2, sd),
                      apply(data, 2, min),
                      apply(data, 2, max))
  colnames(df) <- c('mean', 'sd', 'min', 'max')
  return(print(df))
}

reconstruct_yields <- function(factors, maturities){

  factor_loading_1_vec <- rep(1,length(maturities))
  factor_loading_2_vec <- sapply(maturities/12, function(i) factor_loading_2(i, lambda = 0.0609))
  factor_loading_3_vec <- sapply(maturities/12, function(i) factor_loading_3(i, lambda = 0.0609))
  
  factor_loadings_matrix <- cbind(factor_loading_1_vec, factor_loading_2_vec, factor_loading_3_vec)
  
  reconstructed_yield_row <- function(factor_row){t(factor_loadings_matrix %*% as.numeric(factor_row))}
  reconstructed_yields <- t(apply(factors, 1, reconstructed_yield_row))
                                
  return(reconstructed_yields)
}


load.FRED_MD <- function(){
  
  FRED.MD_data <- read.csv("data/FRED-MD_current.csv", skip = 0, row.names = 1, header =TRUE)[-1,]
  indices <- as.Date(rownames(FRED.MD_data), format = "%m/%d/%Y")
  rownames(FRED.MD_data) <- indices
  
  # columns_na_percentage <- apply(is.na(FRED.MD_data), 2, mean) # percentage of NAs per column
  FRED.MD_data <- apply(FRED.MD_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
  FRED.MD_data <- ts(FRED.MD_data,start=c(1986,1),end=c(2023,12),frequency=12)
  FRED.MD_data <- scale(FRED.MD_data)
  pca_result <- prcomp(FRED.MD_data, center = TRUE, scale. = TRUE)
  macro_factors <- pca_result$x[, 1:3]
  macro_factors <- ts(macro_factors,start=c(1986,1),end=c(2023,12),frequency=12)
  
  return(macro_factors)
}

# load LW data ####
load.LW.adjusted <- function(){
  data = gsheet::gsheet2tbl('https://docs.google.com/spreadsheets/d/1-wmStGZHLx55dSYi3gQK2vb3F8dMw_Nb/edit?usp=drive_link&ouid=117915996921355706819&rtpof=true&sd=true')
  df = data[-(1:8), -1]
  for(i in seq_along(df)) {
    if(is.character(df[[i]])) {
      df[[i]] <- as.numeric(as.character(df[[i]]))
    }
  }
  LW.data = ts(df, start = c(1961,6), frequency=12)
  colnames(LW.data) = 1:360
  LW.data = window(LW.data, start = 1986, end = c(2023,12))
  return(LW.data)
}
