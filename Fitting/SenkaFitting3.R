library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
library(FME)
library(deSolve)


## =======================================================================
## Model growth simulations
## =======================================================================
rm(list = ls())

setwd("/Users/iguex/Documents/CoCulture_Soil/Data/Transfer of data for SynCom model")

Mu_max_vect = matrix(0, nrow = 21, ncol = 96) #Row = species, column = resources.
Mu_max_vect_Monod = matrix(0, nrow = 21, ncol = 96) #Row = species, column = resources.
Lag_time_vect = matrix(0, nrow = 21, ncol = 96) #Row = species, column = resources.
Lag_time_vect_Monod = matrix(0, nrow = 21, ncol = 96) #Row = species, column = resources.
Time_ODE = seq(0,120,1)

for(j in 1:21){
  num_Sheet = j + 6
  Data = read_excel("Senka_Biolog_PlateReader_allstrainscombined.xlsx", sheet = num_Sheet, range = "A1:CU10")
  Names = Data[1,1]
  
  
  Num_iter = length(Data[1,]) - 3
  TT = data.matrix(Data[,2])
  
  for(i in 1:Num_iter){
    
    Res_number = i + 3 
    
    Time_Eval = TT
    
    Data_Cells <- data.matrix(Data[,Res_number]) #Abundance in OD, converted into matrix
    
    diff_OD = abs(Data_Cells[length(Data_Cells)] - Data_Cells[1])
    
    if(diff_OD >= 0.095){
    
      state = c(x = Data_Cells[1])
      
      state_Monod = c(x = Data_Cells[1], R = 3*max(Data_Cells))
      
      x_0 = Data_Cells[1]
      
      Data_Cells <- data.frame(
        time = rep(TT, 1),
        x = c(Data_Cells)
      )
      colnames(Data_Cells) <- c('time','x')
      
      #Function for the logistic estimation
      logist <- function(t, state, parms) {
        with(as.list(c(state, parms)), {
          dx <- 1/(1 + (LT/t)^40)*mu_max*x*(1 - x/Ks) #With lag time
          # dx <- mu_max*x*(1 - x/Ks) #Without lag time
          list(dx)
        })
      }
      
      Monod <- function(t, state, parms){
        with(as.list(c(state, parms)), {
          alpha = 1
          R_conc = max(alpha*R, 0)
          dx = 1/(1 + (LT/t)^40)*x*mu_max*R_conc/(R_conc + Ks)#
          dR = -1/(1 + (LT/t)^40)*x*(mu_max/yield)*R_conc/(R_conc + Ks)
          list(c(dx, dR))
        })
      }
      
      ##===================================
      ## Fitted with logistic model #
      ##===================================
      ## numeric solution 
      ## ODEs system
      parms_init <- c(mu_max = 0.05, Ks = 1.5, LT = 40)
      parms_init_Monod <- c(mu_max = 0.05, Ks = 0.5, LT = 40, yield = 0.3)
      Times <- TT
      
      ## model cost,
      ModelCost2 <- function(P) {
        out <- ode(y = state, func = logist, parms = P, times = TT)
        model = out
        return(modCost(out, Data_Cells)) # object of class modCost
      }
      
      ModelCostMonod <- function(P) {
        out <- ode(y = state_Monod, func = Monod, parms = P, times = TT, atol = 1e-11, rtol = 1e-10)
        model = out
        model = model[,1:2]
        return(modCost(model, Data_Cells)) # object of class modCost
      }
      
      Fit <- modFit(f = ModelCost2, p = parms_init, lower = c(0, x_0, 0),
                    upper = c(2.5, 2, 72))
      
      out <- ode(y = state, func = logist, parms = Fit$par,
                 times = Time_ODE)
      
      
      FitMonod <- modFit(f = ModelCostMonod, p = parms_init_Monod, lower = c(0, 0.01, 0, 0),
                         upper = c(2.5, 2, 72, 1))
      
      outMonod <- ode(y = state_Monod, func = Monod, parms = FitMonod$par,
                  times = Time_ODE, atol = 1e-11, rtol = 1e-10)
      
      Param_est = Fit$par
      Mu_max_vect[j, i] = Param_est[1]
      Lag_time_vect[j, i] = Param_est[3]
      
      Param_est_Monod = FitMonod$par
      
      if(Param_est_Monod[1] > 2.2){
        Mu_max_vect_Monod[j, i] = 0
        Lag_time_vect_Monod[j, i] = 0
      }
      else{
        Mu_max_vect_Monod[j, i] = Param_est_Monod[1]
        Lag_time_vect_Monod[j, i] = Param_est_Monod[3]
      }
  
      pdf(file = paste("/Users/iguex/Documents/CoCulture_Soil/Figures/Fitting figures/Senka_Biolog Monod/", Names, i, ".pdf") ,   # The directory you want to save the file in
          width = 4, # The width of the plot in inches
          height = 4) # The height of the plot in inches
  
      plot(Data_Cells, ylim = c(0, 3), pch = 21, bg = alpha("green", 0.4), col = alpha("green", 0.4))
      lines(out, col = "red", lty = 2)
      lines(outMonod, col = "blue", lty = 2)
      title(paste(Names, i))
      dev.off()
      summary(Fit)
    }
    else{
      x_0 = Data_Cells[1]
      
      Mu_max_vect[j, i] = 0
      Lag_time_vect[j, i] = 0
      
      Mu_max_vect_Monod[j, i] = 0
      Lag_time_vect_Monod[j, i] = 0

      pdf(file = paste("/Users/iguex/Documents/CoCulture_Soil/Figures/Fitting figures/Senka_Biolog Monod/", Names, i, ".pdf") ,   # The directory you want to save the file in
          width = 4, # The width of the plot in inches
          height = 4) # The height of the plot in inches
      
      plot(Data_Cells, ylim = c(0, 3), pch = 21, bg = alpha("green", 0.4), col = alpha("green", 0.4))
      abline(a = x_0, b = 0, col = "red", lty = 2)
      abline(a = x_0, b = 0, col = "blue", lty = 2)
      title(paste(Names, i))
      dev.off()
    }
  }
  Mu_max_vect[j, (Mu_max_vect[j, ] > 3*mean(Mu_max_vect[j, Mu_max_vect[j, ] > 0]))] = 0
  Lag_time_vect_Monod[j, (Mu_max_vect[j, ] > 3*mean(Mu_max_vect[j, Mu_max_vect[j, ] > 0]))] = 0
}

Mu_max_vect = as.matrix(Mu_max_vect)
clipr::write_clip(Mu_max_vect)

Lag_time_vect = as.matrix(Lag_time_vect)
clipr::write_clip(Lag_time_vect)

Mu_max_vect_Monod = as.matrix(Mu_max_vect_Monod)
clipr::write_clip(Mu_max_vect_Monod)

Lag_time_vect_Monod = as.matrix(Lag_time_vect_Monod)
clipr::write_clip(Lag_time_vect_Monod)