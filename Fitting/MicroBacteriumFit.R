library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
library(Hmisc)
library(FME)
library(deSolve)
library(MASS)


## =======================================================================
## Model growth simulations
## =======================================================================
rm(list = ls())

setwd("/Users/isalinelucille-guex/Documents/CoCulture_Soil/Data/Transfer of data for SynCom model")

Data = read_excel("microbacteriumgrowthcurve.xlsx", sheet = 1, range = "A1:F6") #Lag time Monod
Mu_max_vect_Monod = matrix(0, nrow = 1, ncol = 3) #Row = species, column = resources.
Lag_time_vect_Monod = matrix(0, nrow = 1, ncol = 3) #Row = species, column = resources.
Time_ODE = seq(0,74,1)

Num_iter = 3
TT = data.matrix(Data[,1])
Names = "microbacterium"  
for(i in 1:Num_iter){
    
  Time_Eval = TT
    
  Data_Cells <- exp(data.matrix(Data[,i+1])) #Abundance in OD, converted into matrix
  state_Monod = c(x = Data_Cells[1], R = 3*max(Data_Cells))
  x_0 = Data_Cells[1]
  Data_Cells <- data.frame(
    time = rep(TT, 1),
    x = c(Data_Cells)
    )
  colnames(Data_Cells) <- c('time','x')
      
  Monod <- function(t, state, parms){
    with(as.list(c(state, parms)), {
      alpha = 1
      R_conc = max(alpha*R, 0)
      dx = 1/(1 + (LT/t)^40)*x*mu_max*R_conc/(R_conc + Ks)#
      dR = -1/(1 + (LT/t)^40)*x*(mu_max/yield)*R_conc/(R_conc + Ks)
      list(c(dx, dR))
     })
  }
      
  parms_init_Monod <- c(mu_max = 0.05, Ks = 50, LT = 4, yield = 0.3)
  Times <- TT
      
  ModelCostMonod <- function(P) {
    out <- ode(y = state_Monod, func = Monod, parms = P, times = TT, atol = 1e-11, rtol = 1e-10)
    model = out
    model = model[,1:2]
    return(modCost(model, Data_Cells)) # object of class modCost
  }
      
      
  FitMonod <- modFit(f = ModelCostMonod, p = parms_init_Monod, lower = c(0, 0.01, 0, 0),
                      upper = c(2.5, 5000, 72, 1))
      
  outMonod <- ode(y = state_Monod, func = Monod, parms = FitMonod$par,
                  times = Time_ODE, atol = 1e-11, rtol = 1e-10)
      
  Param_est_Monod = FitMonod$par
  Mu_max_vect_Monod[i] = Param_est_Monod[1]
  Lag_time_vect_Monod[i] = Param_est_Monod[3]
      
      
  plot(Data_Cells, xlim = c(0, max(TT)), pch = 21, bg = alpha("green", 0.4), col = alpha("green", 0.4))
  lines(outMonod, col = "blue", lty = 2)
  title(paste(Names, i))
  summary(FitMonod)
}

Mu_max_vect = as.matrix(Mu_max_vect)
clipr::write_clip(Mu_max_vect)

Lag_time_vect = as.matrix(Lag_time_vect)
clipr::write_clip(Lag_time_vect)

Mu_max_vect_Monod = as.matrix(Mu_max_vect_Monod)
clipr::write_clip(Mu_max_vect_Monod)

Lag_time_vect_Monod = as.matrix(Lag_time_vect_Monod)
clipr::write_clip(Lag_time_vect_Monod)