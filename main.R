# Load packages
source('R/setup.R')

source('R/utility.R')
source('R/fcst.AFFM.R')
source('R/fcst.DNS.R')
source('R/load.FRED-MD.R')
source('R/generate.PC.factors.R')
source('R/oos.results.R')
source('R/oos.results.RW.R')

### Preprocessing etc. ####

# Load data
yields = load.LW.adjusted() # load data with adjusted function to fix old link
FRED.MD_data <- load.FRED_MD()

# Start Process
workinggrid = seq(from=as.numeric(colnames(yields)[1]),to=as.numeric(rev(colnames(yields))[1]), by=0.5) # Define working grid
fdaobj = fda.preprocess(yields, workinggrid = workinggrid)
observationgrid = fdaobj$observationgrid

yield_data <- vector("list", max(456))
IC <- vector("list", max(456))
macro_factors_3 <- vector("list", max(456))
macro_factors_2 <- vector("list", max(456))
macro_factors_1 <- vector("list", max(456))
for (i in 120:456) {
  yields_window <- window(yields, start = time(yields)[1], end = time(yields)[i])
  fdaobj_current = fda.preprocess(yields_window, workinggrid = workinggrid)
  yield_data[[i]] = fts.cumAC(fdaobj_current)
  IC[[i]] <- fts.criterion(fdaobj_current, K.max = 12, p.max = 12)$IC.min
  FRED.MD_window <- window(FRED.MD_data, start = time(yields)[1], end = time(yields)[i])
  macro_factors_3[[i]] <- generate.PC.factors(data = FRED.MD_window, no.factors = 3)
  macro_factors_2[[i]] <- generate.PC.factors(data = FRED.MD_window, no.factors = 2)
  macro_factors_1[[i]] <- generate.PC.factors(data = FRED.MD_window, no.factors = 1)
  cat("\r",i)
}
 save(yield_data, file = "results/yield_data.RData")
 save(IC, file = "results/IC.RData")
 save(macro_factors_3, file = "results/macro_factors_3.RData")
 save(macro_factors_2, file = "results/macro_factors_2.RData")
 save(macro_factors_1, file = "results/macro_factors_1.RData")

load("results/yield_data.RData")
load("results/IC.RData")
load("results/macro_factors_3.RData")
load("results/macro_factors_2.RData")
load("results/macro_factors_1.RData")

### Run Models ####

model_AR.3 <- list(func = fcst.AFFM, args = list(K = 3, p = 1, AR = TRUE))
model_AR.4 <- list(func = fcst.AFFM, args = list(K = 4, p = 1, AR = TRUE))
model_AR.6 <- list(func = fcst.AFFM, args = list(K = 6, p = 1, AR = TRUE))
model_AR.8 <- list(func = fcst.AFFM, args = list(K = 8, p = 1, AR = TRUE))
model_VAR <- list(func = fcst.AFFM, args = list(AR = FALSE))

model_DNS_AR <- list(func = fcst.DNS, args = list(p = 1, AR = TRUE))
model_DNS_VAR <- list(func = fcst.DNS, args = list(p = 1, AR = FALSE))

# Baseline
results_AR.3 <- oos.results(model = model_AR.3, yield_data = yield_data)
  #save(results_AR.3, file = "results/results_AR.3.RData")
  #load("results/results_AR.3.RData")
results_AR.4 <- oos.results(model = model_AR.4, yield_data = yield_data)
  #save(results_AR.4, file = "results/results_AR.4.RData")
  #load("results/results_AR.4.RData")
results_AR.6 <- oos.results(model = model_AR.6, yield_data = yield_data)
  #save(results_AR.6, file = "results/results_AR.6.RData")
  #load("results/results_AR.6.RData")
results_AR.8 <- oos.results(model = model_AR.8, yield_data = yield_data)
  #save(results_AR.8, file = "results/results_AR.8.RData")
  #load("results/results_AR.8.RData")

results_AR.BIC <- oos.results(model = model_AR.3, yield_data = yield_data, IC = IC, IC_type = "BIC")
  #save(results_AR.BIC, file = "results/results_AR.BIC.RData")
  #load("results/results_AR.BIC.RData")
results_AR.HQC <- oos.results(model = model_AR.3, yield_data = yield_data, IC = IC, IC_type = "HQC")
  #save(results_AR.HQC, file = "results/results_AR.HQC.RData")
  #load("results/results_AR.HQC.RData")
results_VAR.BIC <- oos.results(model = model_VAR, yield_data = yield_data, IC = IC, IC_type = "BIC")
  #save(results_VAR.BIC, file = "results/results_VAR.BIC.RData")
  #load("results/results_VAR.BIC.RData")
results_VAR.HQC <- oos.results(model = model_VAR, yield_data = yield_data, IC = IC, IC_type = "HQC")
  #save(results_VAR.HQC, file = "results/results_VAR.HQC.RData")
  #load("results/results_VAR.HQC.RData")

#### Macro PC3 ####
results_AR.3.PC3 <- oos.results(model = model_AR.3, yield_data = yield_data, macro_factors = macro_factors_3)
  #save(results_AR.3.PC3, file = "results/results_AR.3.PC3.RData")
  #load("results/results_AR.3.PC3.RData")
results_AR.4.PC3 <- oos.results(model = model_AR.4, yield_data = yield_data, macro_factors = macro_factors_3)
  #save(results_AR.4.PC3, file = "results/results_AR.4.PC3.RData")
  #load("results/results_AR.4.PC3.RData")
results_AR.6.PC3 <- oos.results(model = model_AR.6, yield_data = yield_data, macro_factors = macro_factors_3)
  #save(results_AR.6.PC3, file = "results/results_AR.6.PC3.RData")
  #load("results/results_AR.6.PC3.RData")
results_AR.8.PC3 <- oos.results(model = model_AR.8, yield_data = yield_data, macro_factors = macro_factors_3)
  #save(results_AR.8.PC3, file = "results/results_AR.8.PC3.RData")
  #load("results/results_AR.8.PC3.RData")

results_AR.BIC.PC3 <- oos.results(model = model_AR.3, yield_data = yield_data, IC = IC, IC_type = "BIC", macro_factors = macro_factors_3)
  #save(results_AR.BIC.PC3, file = "results/results_AR.BIC.PC3.RData")
  #load("results/results_AR.BIC.PC3.RData")
results_AR.HQC.PC3 <- oos.results(model = model_AR.3, yield_data = yield_data, IC = IC, IC_type = "HQC", macro_factors = macro_factors_3)
  #save(results_AR.HQC.PC3, file = "results/results_AR.HQC.PC3.RData")
  #load("results/results_AR.HQC.PC3.RData")
results_VAR.BIC.PC3 <- oos.results(model = model_VAR, yield_data = yield_data, IC = IC, IC_type = "BIC", macro_factors = macro_factors_3)
  #save(results_VAR.BIC.PC3, file = "results/results_VAR.BIC.PC3.RData")
  #load("results/results_VAR.BIC.PC3.RData")
results_VAR.HQC.PC3 <- oos.results(model = model_VAR, yield_data = yield_data, IC = IC, IC_type = "HQC", macro_factors = macro_factors_3)
  #save(results_VAR.HQC.PC3, file = "results/results_VAR.HQC.PC3.RData")
  #load("results/results_VAR.HQC.PC3.RData")

##### Macro PC2 ####
results_AR.3.PC2 <- oos.results(model = model_AR.3, yield_data = yield_data, macro_factors = macro_factors_2)
  #save(results_AR.3.PC2, file = "results/results_AR.3.PC2.RData")
  #load("results/results_AR.3.PC2.RData")
results_AR.4.PC2 <- oos.results(model = model_AR.4, yield_data = yield_data, macro_factors = macro_factors_2)
  #save(results_AR.4.PC2, file = "results/results_AR.4.PC2.RData")
  #load("results/results_AR.4.PC2.RData")
results_AR.6.PC2 <- oos.results(model = model_AR.6, yield_data = yield_data, macro_factors = macro_factors_2)
  #save(results_AR.6.PC2, file = "results/results_AR.6.PC2.RData")
  #load("results/results_AR.6.PC2.RData")
results_AR.8.PC2 <- oos.results(model = model_AR.8, yield_data = yield_data, macro_factors = macro_factors_2)
  #save(results_AR.8.PC2, file = "results/results_AR.8.PC2.RData")
  #load("results/results_AR.8.PC2.RData")

results_AR.BIC.PC2 <- oos.results(model = model_AR.3, yield_data = yield_data, IC = IC, IC_type = "BIC", macro_factors = macro_factors_2)
  #save(results_AR.BIC.PC2, file = "results/results_AR.BIC.PC2.RData")
  #load("results/results_AR.BIC.PC2.RData")
results_AR.HQC.PC2 <- oos.results(model = model_AR.3, yield_data = yield_data, IC = IC, IC_type = "HQC", macro_factors = macro_factors_2)
  #save(results_AR.HQC.PC2, file = "results/results_AR.HQC.PC2.RData")
  #load("results/results_AR.HQC.PC2.RData")
results_VAR.BIC.PC2 <- oos.results(model = model_VAR, yield_data = yield_data, IC = IC, IC_type = "BIC", macro_factors = macro_factors_2)
  #save(results_VAR.BIC.PC2, file = "results/results_VAR.BIC.PC2.RData")
  #load("results/results_VAR.BIC.PC2.RData")
results_VAR.HQC.PC2 <- oos.results(model = model_VAR, yield_data = yield_data, IC = IC, IC_type = "HQC", macro_factors = macro_factors_2)
  #save(results_VAR.HQC.PC2, file = "results/results_VAR.HQC.PC2.RData")
  #load("results/results_VAR.HQC.PC2.RData")

#### Macro PC1 ####
results_AR.3.PC1 <- oos.results(model = model_AR.3, yield_data = yield_data, macro_factors = macro_factors_1)
  #save(results_AR.3.PC1, file = "results/results_AR.3.PC1.RData")
  #load("results/results_AR.3.PC1.RData")
results_AR.4.PC1 <- oos.results(model = model_AR.4, yield_data = yield_data, macro_factors = macro_factors_1)
  #save(results_AR.4.PC1, file = "results/results_AR.4.PC1.RData")
  #load("results/results_AR.4.PC1.RData")
results_AR.6.PC1 <- oos.results(model = model_AR.6, yield_data = yield_data, macro_factors = macro_factors_1)
  #save(results_AR.6.PC1, file = "results/results_AR.6.PC1.RData")
  #load("results/results_AR.6.PC1.RData")
results_AR.8.PC1 <- oos.results(model = model_AR.8, yield_data = yield_data, macro_factors = macro_factors_1)
  #save(results_AR.8.PC1, file = "results/results_AR.8.PC1.RData")
  #load("results/results_AR.8.PC1.RData")

results_AR.BIC.PC1 <- oos.results(model = model_AR.3, yield_data = yield_data, IC = IC, IC_type = "BIC", macro_factors = macro_factors_1)
  #save(results_AR.BIC.PC1, file = "results/results_AR.BIC.PC1.RData")
  #load("results/results_AR.BIC.PC1.RData")
results_AR.HQC.PC1 <- oos.results(model = model_AR.3, yield_data = yield_data, IC = IC, IC_type = "HQC", macro_factors = macro_factors_1)
  #save(results_AR.HQC.PC1, file = "results/results_AR.HQC.PC1.RData")
  #load("results/results_AR.HQC.PC1.RData")
results_VAR.BIC.PC1 <- oos.results(model = model_VAR, yield_data = yield_data, IC = IC, IC_type = "BIC", macro_factors = macro_factors_1)
  #save(results_VAR.BIC.PC1, file = "results/results_VAR.BIC.PC1.RData")
  #load("results/results_VAR.BIC.PC1.RData")
results_VAR.HQC.PC1 <- oos.results(model = model_VAR, yield_data = yield_data, IC = IC, IC_type = "HQC", macro_factors = macro_factors_1)
  #save(results_VAR.HQC.PC1, file = "results/results_VAR.HQC.PC1.RData")
  #load("results/results_VAR.HQC.PC1.RData")


# DNS
results_DNS_AR <- oos.results(model = model_DNS_AR, yield_data = yield_data)
  #save(results_DNS_AR, file = "results/results_DNS_AR.RData")
  #load("results/results_DNS_AR.RData")
results_DNS_VAR <- oos.results(model = model_DNS_VAR, yield_data = yield_data)
  #save(results_DNS_VAR, file = "results/results_DNS_VAR.RData")
  #load("results/results_DNS_VAR.RData")

results_DNS_AR.PC1 <- oos.results(model = model_DNS_AR, yield_data = yield_data, macro_factors = macro_factors_1)
  #save(results_DNS_AR.PC1, file = "results/results_DNS_AR.PC1.RData")
  #load("results/results_DNS_AR.PC1.RData")
results_DNS_VAR.PC1 <- oos.results(model = model_DNS_VAR, yield_data = yield_data, macro_factors = macro_factors_1)
  #save(results_DNS_VAR.PC1, file = "results/results_DNS_VAR.PC1.RData")
  #load("results/results_DNS_VAR.PC1.RData")

results_DNS_AR.PC2 <- oos.results(model = model_DNS_AR, yield_data = yield_data, macro_factors = macro_factors_2)
  #save(results_DNS_AR.PC2, file = "results/results_DNS_AR.PC2.RData")
  #load("results/results_DNS_AR.PC2.RData")
results_DNS_VAR.PC2 <- oos.results(model = model_DNS_VAR, yield_data = yield_data, macro_factors = macro_factors_2)
  #save(results_DNS_VAR.PC2, file = "results/results_DNS_VAR.PC2.RData")
  #load("results/results_DNS_VAR.PC2.RData")

results_DNS_AR.PC3 <- oos.results(model = model_DNS_AR, yield_data = yield_data, macro_factors = macro_factors_3)
  #save(results_DNS_AR.PC3, file = "results/results_DNS_AR.PC3.RData")
  #load("results/results_DNS_AR.PC3.RData")
results_DNS_VAR.PC3 <- oos.results(model = model_DNS_VAR, yield_data = yield_data, macro_factors = macro_factors_3)
  #save(results_DNS_VAR.PC3, file = "results/results_DNS_VAR.PC3.RData")
  #load("results/results_DNS_VAR.PC3.RData")

results_RW <- oos.results.RW(yield_data = yield_data)


# Show results of AR, K = 3, p = 1 as example
horizons = c(1,3,6,12)
maturities <- c(1,3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120,180,240,300,360)

results_AR.3$output.overall

model <- results_AR.3
output_RMSFE_overall <- round(c(model$insample.fit.avg[2], as.vector(model$output.overall[,1])),3)
output_RMSFE_maturities <- round(do.call(rbind, lapply(model$RMSFE.maturities[horizons], `[`, maturities)), 3)
output_coverage_overall <- round(as.vector(model$output.overall[,2]),3)
output_coverage_overall_GARCH <- round(as.vector(model$output.overall[,3]),3)
output_coverage_maturities <- round(do.call(rbind, lapply(model$coverage.maturities[horizons], `[`, match(observationgrid, workinggrid)[maturities])), 3)
output_coverage_maturities_GARCH <- round(do.call(rbind, lapply(model$coverage.GARCH.maturities[horizons], `[`, match(observationgrid, workinggrid)[maturities])), 3)


### Generate Figures
source("scripts/figures.R")


