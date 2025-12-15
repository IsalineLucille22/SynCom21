library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
# library(Hmisc)
library(FME)
library(deSolve)
library(MASS)


## =======================================================================
## Model growth simulations
## =======================================================================
rm(list = ls())

#setwd("/Users/isalinelucille-guex/Documents/CoCulture_Soil/Data/Transfer of data for SynCom model")
setwd("/Users/iguex/Documents/CoCulture_Soil/Data/Transfer of data for SynCom model")
# Data = read_excel("PTYG_OD600_21strain_Senka&Tania.xlsx", sheet = 5, range = "A74:DM94") #Mu_max without lag time
# Data = read_excel("PTYG_OD600_21strain_Senka&Tania.xlsx", sheet = 5, range = "A192:DI212") #Mu_max with lag time
# Data = read_excel("PTYG_OD600_21strain_Senka&Tania.xlsx", sheet = 5, range = "A238:DI258") #Lag time logistic
Data = read_excel("PTYG_OD600_21strain_Senka&Tania.xlsx", sheet = 7, range = "A213:DI233") #Mu_max Monod
# Data = read_excel("PTYG_OD600_21strain_Senka&Tania.xlsx", sheet = 7, range = "A259:DI279") #Lag time Monod
# Data = read_excel("Senka_Biolog_PlateReader_allstrainscombined.xlsx", sheet = 1, range = "A213:DI233") 

# #With Lysobacter
# Data = read_excel("PTYG_OD600_21strain_Senka&Tania.xlsx", sheet = 9, range = "A75:CV96") #CV for Microbacterium only with sheet = 9, DI for the other Mu_max Monod with sheet = 8
# Data = read_excel("PTYG_OD600_21strain_Senka&Tania.xlsx", sheet = 9, range = "A124:CV145") #CV for Microbacterium only, DI for the other Lag time Monod


Nb_species = length(data.matrix(Data[,2]))
Nb_res = length(Data[1,]) - 1

yield_vect = matrix(0, nrow = Nb_species, ncol = 1) 

parameters_norm = matrix(0, nrow = Nb_species, ncol = 2) 

for(i in 1:Nb_species){
  # i = 16
  Growth_rates_tot = data.matrix(Data[i,2:100])#Data[i,2:116] without lag time, Data[i,2:113] with lag time. Data[i,2:100] for Microbacterium only, Data[i,2:113] for the other species
  na_val = is.na(Growth_rates_tot)
  Growth_rates_tot = Growth_rates_tot[!na_val] #Remove
  Growth_rates = Growth_rates_tot[Growth_rates_tot > 0]
  yield = length(Growth_rates)/length(Growth_rates_tot)
  yield_vect[i] = yield
  # fit <- fitdistr(Growth_rates, "normal")
  # fit <- fitdistr(Growth_rates, "weibull")
  # fit <- fitdistr(Growth_rates, "gamma")
  fit <- fitdistr(Growth_rates, "log-normal")
  para <- fit$estimate
  parameters_norm[i,] = para
  x2 <- seq(0, 2*max(Growth_rates), length = 1400)
  # Log-Normal curve
  fun <- dlnorm(x2, meanlog = para[1], sdlog = para[2])
  # pdf(file = paste("/Users/iguex/Documents/CoCulture_Soil/Figures/Fitting figures/Fit_LT/", Data[i,1], ".pdf") ,   # The directory you want to save the file in
  #     width = 4, # The width of the plot in inches
  #     height = 4) # The height of the plot in inches
  hist(Growth_rates, main = Data[i,1], prob = TRUE, ylim = c(0, 1.2*max(fun)), xlim = c(0, 1.1*max(Growth_rates)))
  lines(x2, fun, col = 2, lwd = 2)
  # curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
  # curve(dweibull(x, shape = para[1], scale = para[2]), col = 2, add = TRUE)
  # curve(dgamma(x, shape = para[1], rate = para[2]), col = 2, add = TRUE)
  # curve(dlnorm(x, meanlog = para[1], sdlog = para[2]), col = 2, add = TRUE)
  title(Data[i,1])
  # dev.off()
}

parameters_norm = as.matrix(parameters_norm)
clipr::write_clip(parameters_norm)

na_val_vect = as.matrix(yield_vect)
clipr::write_clip(yield_vect)

