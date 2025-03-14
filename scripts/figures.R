# Load packages
source('R/setup.R')

### Generate figures ####
yields = load.LW.adjusted() # load data with adjusted function to fix old link
FRED.MD_data <- load.FRED_MD()
load("results/yield_data.RData")

workinggrid = seq(from=as.numeric(colnames(yields)[1]),to=as.numeric(rev(colnames(yields))[1]), by=0.5) # Define working grid
yields_fct <- fdaobj$densedata

library(plotly)
library(magrittr)
library(ggplot2)
library(reshape2)

time <- seq(1986, 2023, length.out = nrow(yields_fct))  
maturities <- workinggrid                

specific_ticks <- c(1, 60, 120, 180, 240, 300) 
specific_labels <- as.character(specific_ticks)        

plot_ly(
  x = maturities,     
  y = time,           
  z = as.matrix(yields_fct), 
  type = "surface",   
  colors = "YlGnBu",  
  colorbar = list(
    title = "Yield (%)",  
    len = 0.5,            
    thickness = 15        
  ) 
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = "Maturity (Months)",
        tickmode = "array",          
        tickvals = specific_ticks,   
        ticktext = specific_labels   
      ),
      yaxis = list(title = "Time"),
      zaxis = list(title = "Yield (%)"),
      eye = list(x = 1.5, y = 1.5, z = 0.8),
      aspectmode = "manual",          
      aspectratio = list(x = 1, y = 1.5, z = 1) 
    )
  ) # Has to be saved manually !


# DNS loadings
slope <- function(m, lambda = 0.0609){(1-exp(-lambda*m))/(lambda * m)}
curvature <- function(m, lambda = 0.0609){(1-exp(-lambda*m))/(lambda *m) - exp(- lambda * m)}

loadings_DNS <- t(sapply(workinggrid, function(i) {
  c(1, slope(i), curvature(i))
}))
loadings_DNS <- as.data.frame(loadings_DNS)
colnames(loadings_DNS) <- c("Level", "Slope", "Curvature")
data_DNS <- melt(loadings_DNS, variable.name = "Factor", value.name = "Value")
data_DNS$Maturity <- rep(workinggrid, times = 3)

ggplot(data_DNS, aes(x = Maturity, y = Value, color = Factor, linetype = Factor)) +
  geom_line(size = 1.1) +  
  scale_color_manual(
    values = c("black", "grey20", "grey40", "grey60"),
    labels = c(expression(level), expression(slope), expression(curvature))# Latex-Labels
  ) +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted"),
    labels = c(expression(level), expression(slope), expression(curvature)) # Latex-Labels
  ) +
  scale_y_continuous(limits = c(0, 1.2)) +
  labs(
    title = "DNS loading functions",
    x = "Maturity (months)",
    y = "",
    color = "Functions", 
    linetype = "Functions"
  ) +
  theme_minimal(base_size = 14) +  
  theme(
    panel.grid = element_blank(),                                                
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    plot.title = element_text(hjust = 0.5, face = "bold"),    
    plot.subtitle = element_text(hjust = 0.5, face = "italic"),
    legend.position = "top",                                  
    legend.title = element_blank(),                           
    legend.text = element_text(size = 12),                    
    legend.key.width = unit(2, "cm"),                         
    legend.spacing.x = unit(0.5, "cm")                    
  )
ggsave(
  filename = "figures/loadings_DNS.png",  
  plot = last_plot(),                
  width = 12,                        
  height = 6,                      
  dpi = 300                         
)


# AFFM loadings T = 120
loadings_AFFM_120 <- yield_data[[120]]$eigenfunctions[,1:4,drop=F]
data_AFFM_120 <- melt(loadings_AFFM_120, value.name = "Value")
data_AFFM_120$Maturity <- rep(workinggrid, times = 4)

ggplot(data_AFFM_120, aes(x = Maturity, y = Value, color = Var2, linetype = Var2)) +
  geom_line(size = 1.1) +  
  scale_color_manual(
    values = c("black", "grey20", "grey40", "grey60"),
    labels = c(expression(delta[1]), expression(delta[2]), expression(delta[3]), expression(delta[4]))# Latex-Labels
  ) +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted", "dotdash"),
    labels = c(expression(delta[1]), expression(delta[2]), expression(delta[3]), expression(delta[4])) # Latex-Labels
  ) +
  labs(
    title = "First four AFFM loading functions, based on sample size T = 120",
    subtitle = "(January 1986 - December 1995)",
    x = "Maturity (months)",
    y = "",
    color = "Functions",  # Legendentitel
    linetype = "Functions"
  ) +
  theme_minimal(base_size = 14) +  
  theme(
    panel.grid = element_blank(),
    #panel.grid.major = element_line(color = "grey80", size = 0.5), 
    #panel.grid.minor = element_blank(),                                                  
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    plot.title = element_text(hjust = 0.5, face = "bold"),    
    plot.subtitle = element_text(hjust = 0.5, face = "italic"), 
    legend.position = "top",                                  
    legend.title = element_blank(),                           
    legend.text = element_text(size = 12),                    
    legend.key.width = unit(2, "cm"),                         
    legend.spacing.x = unit(0.5, "cm")                    
  )
ggsave(
  filename = "figures/loadings_AFFM_120.png",  
  plot = last_plot(),                
  width = 12,                        
  height = 6,                      
  dpi = 300                         
)

# AFFM loadings T = 456
loadings_AFFM_456 <- yield_data[[456]]$eigenfunctions[,1:4,drop=F]
data_AFFM_456 <- melt(loadings_AFFM_456, value.name = "Value")
data_AFFM_456$Maturity <- rep(workinggrid, times = 4)

ggplot(data_AFFM_456, aes(x = Maturity, y = Value, color = Var2, linetype = Var2)) +
  geom_line(size = 1.1) +  
  scale_color_manual(
    values = c("black", "grey20", "grey40", "grey60"),
    labels = c(expression(delta[1]), expression(delta[2]), expression(delta[3]), expression(delta[4]))# Latex-Labels
  ) +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted", "dotdash"),
    labels = c(expression(delta[1]), expression(delta[2]), expression(delta[3]), expression(delta[4])) # Latex-Labels
  ) +
  labs(
    title = "First four AFFM loading functions, based on sample size T = 456",
    subtitle = "(January 1986 - December 2023)",
    x = "Maturity (months)",
    y = "",
    color = "Functions",  # Legendentitel
    linetype = "Functions"
  ) +
  theme_minimal(base_size = 14) +  # Minimalistisches Theme mit größerer Schrift
  theme(
    panel.grid = element_blank(),
    #panel.grid.major = element_line(color = "grey80", size = 0.5), 
    #panel.grid.minor = element_blank(),                                                   
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    plot.title = element_text(hjust = 0.5, face = "bold"),    
    plot.subtitle = element_text(hjust = 0.5, face = "italic"), 
    legend.position = "top",                                  
    legend.title = element_blank(),                           
    legend.text = element_text(size = 12),                    
    legend.key.width = unit(2, "cm"),                         
    legend.spacing.x = unit(0.5, "cm")                    
  )
ggsave(
  filename = "figures/loadings_AFFM_456.png",  
  plot = last_plot(),                
  width = 12,                        
  height = 6,                      
  dpi = 300                         
)