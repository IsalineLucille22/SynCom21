# install.packages("Hmisc", repos = "https://cloud.r-project.org/")
library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
library(FME)
library(deSolve)
#library(Hmisc)
# install.packages(c("FME", "deSolve", "rootSolve", "coda"))





## =======================================================================
## Model growth simulations
## =======================================================================
rm(list = ls())

setwd("/Users/iguex/Documents/CoCulture_Soil/Data/Transfer of data for SynCom model")

Data = read_excel("PTYG_OD600_21strain_Senka&Tania.xlsx", sheet = 3, range = "A1:CG233")
Names = read_excel("PTYG_OD600_21strain_Senka&Tania.xlsx", sheet = 4, range = "A1:D85")

Num_iter = length(Data[1,]) - 1
TT = data.matrix(Data[,1])
Mu_max_vect = c()
Mu_max_vect_Monod = c()
Mu_max_vect_Monod_SANN = c()
Lag_time_vect = c()
Lag_time_vect_Monod = c()
Lag_time_vect_Monod_SANN = c()
nb_obs = length(data.matrix(Data[,1]))


for(i in 1:Num_iter){#Num_iter
  Species_number = i + 1
  
  if(Data[1, Species_number] > Data[nb_obs, Species_number]){
    Data_temp = Data[1:150,]
  }
  else{
    Data_temp = Data
  }
    
  nb_data_pts = length(data.matrix(Data_temp[,1]))
  TT = data.matrix(Data_temp[,1])
  Time_Eval = TT
  
  Data_Cells <- data.matrix(Data_temp[,Species_number]) #Abundance in OD, converted into matrix
  
  x_0 = Data_Cells[1]
  
  diff_OD = abs(max(Data_Cells) - min(Data_Cells))
  
  if(diff_OD >= 0.02){
  
    state = c(x = mean(Data_Cells[1,]))
    
    R = min(3*max(Data_Cells), 3*mean(max(Data_Cells), max(Data_Cells[nb_data_pts],0)))
    state_Monod = c(x = mean(Data_Cells[1,]), R = R)
    
    Data_Cells <- data.frame(
      time = TT,
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
        dx = 1/(1 + (LT/t)^40)*x*mu_max*R_conc/(R_conc + Ks)
        dR = -1/(1 + (LT/t)^40)*x*(mu_max/yield)*R_conc/(R_conc + Ks)
        list(c(dx, dR))
      })
    }
    
    ##===================================
    ## Fitted with logistic model #
    ##===================================
    ## numeric solution 
    ## ODEs system
    parms_init <- c(mu_max = 0.5, Ks = 0.1, LT = 40)
    parms_init_Monod <- c(mu_max = 0.5, Ks = 0.5, LT = 40, yield = 0.3)
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
    
    
    Fit <- modFit(f = ModelCost2, p = parms_init, lower = c(0, 0, 0),
                  upper = c(2.5, 2, 72))
    
    out <- ode(y = state, func = logist, parms = Fit$par,
               times = Time_Eval)
    
    FitMonod_SANN <- modFit(f = ModelCostMonod, p = parms_init_Monod, method = "SANN", lower = c(0, 0, 0, 0),
                       upper = c(2.5, 2, 72, 1))
    
    FitMonod <- modFit(f = ModelCostMonod, p = parms_init_Monod, method = "SANN", lower = c(0, 0, 0, 0),
                       upper = c(2.5, 2, 72, 1))
    
    outMonod_SANN <- ode(y = state_Monod, func = Monod, parms = FitMonod_SANN$par,
                    times = Time_Eval, atol = 1e-11, rtol = 1e-10)
    
    outMonod <- ode(y = state_Monod, func = Monod, parms = FitMonod$par,
                    times = Time_Eval, atol = 1e-11, rtol = 1e-10)
    
    Param_est = Fit$par
    Mu_max_vect[i] = Param_est[1]
    Lag_time_vect[i] = Param_est[3]
    
    Param_est_Monod = FitMonod$par
    Mu_max_vect_Monod[i] = Param_est_Monod[1]
    Lag_time_vect_Monod[i] = Param_est_Monod[3]
    
    Param_est_Monod_SANN = FitMonod_SANN$par
    Mu_max_vect_Monod_SANN[i] = Param_est_Monod_SANN[1]
    Lag_time_vect_Monod_SANN[i] = Param_est_Monod_SANN[3]
    
    pdf(file = paste("/Users/iguex/Documents/CoCulture_Soil/Figures/Fitting figures/Tania&Senka Monod/", Names[i,1], i, "test.pdf") ,   # The directory you want to save the file in
        width = 4, # The width of the plot in inches
        height = 4) # The height of the plot in inches
  
    plot(Data_Cells, xlim = c(0, max(TT)), pch = 21, bg = alpha("green", 0.4), col = alpha("green", 0.4))
    lines(out, col = "red", lty = 2)
    lines(outMonod, col = "blue", lty = 2)
    lines(outMonod_SANN, col = "pink", lty = 2)
    title(Names[i,1])
    dev.off()
    summary(FitMonod)
  }
  else{
    Mu_max_vect[i] = 0
    Lag_time_vect[i] = 0
    
    Mu_max_vect_Monod[i] = 0 
    Lag_time_vect_Monod[i] = 0

    pdf(file = paste("/Users/iguex/Documents/CoCulture_Soil/Figures/Fitting figures/Tania&Senka Monod/", Names[i,1], i, "test.pdf") ,   # The directory you want to save the file in
        width = 4, # The width of the plot in inches
        height = 4) # The height of the plot in inches

    plot(Data_Cells, xlim = c(0, max(TT)), pch = 21, bg = alpha("green", 0.4), col = alpha("green", 0.4))
    abline(a = x_0, b = 0, col = "red", lty = 2)
    abline(a = x_0, b = 0, col = "blue", lty = 2)
    title(Names[i,1])
    dev.off()
  }
}

Mu_max_vect = as.matrix(Mu_max_vect)
clipr::write_clip(Mu_max_vect)

Lag_time_vect = as.matrix(Lag_time_vect)
clipr::write_clip(Lag_time_vect)

Mu_max_vect_Monod = as.matrix(Mu_max_vect_Monod)
clipr::write_clip(Mu_max_vect_Monod)

Lag_time_vect_Monod = as.matrix(Lag_time_vect_Monod)
clipr::write_clip(Lag_time_vect_Monod)
