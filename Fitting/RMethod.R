library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
library(Hmisc)
library(FME)
library(deSolve)


## =======================================================================
## Model growth simulations
## =======================================================================
rm(list = ls())

setwd("/Users/isalinelucille-guex/Documents/CoCulture_Soil/Data")

Data = read_excel("Phil_S20_community_growth.xlsx", sheet = 4, range = "A24:G44")

Species_number = 10 #Values from 1 to 20

#TT <- c(0,0,0,0,1,1,1,1,3,3,3,3,7,7,7,7,15,15,15,15)*24#,22,22,22,22)*24 #Time converted into hours
TT <- c(0,1,3,7,15)*24 #Time converted into hours


Time_Eval = seq(0, 480, 1)

#Data_Cells <- Data[Species_number,c(2:21)]
Data_Cells <- Data[Species_number,c(3:7)]

Data_Cells = as.numeric(Data_Cells)
Data_Cells = Data_Cells*1e15
R_0 = 1.5*1e-03*1e15
state = c(x = Data_Cells[1])
state_1 = c(x = Data_Cells[1], y = 0, c = 0, r = R_0, z = Data_Cells[1])
#state_sts = c(x = Data_Cells[1])
state_sts = c(x = Data_Cells[1], y = 0, c = 0, r = R_0, z = Data_Cells[1])

Data_Cells <- data.frame(
  time = TT,
  x = Data_Cells
)

#Function for the logistic estimation
logist <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dx <- mu_max*x*(1 - x/(R_0*rho))
    list(dx)
  })
}

#Function for Monod estimation
Monod <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    dx = -kappa_1*x*r + (2*kappa_2 + kappa_3)*y
    dy = kappa_1*x*r - (kappa_2 + kappa_3)*y
    dc = kappa_3*y
    dr = -kappa_1*x*r
    dz = dx + dy
    list(c(dx, dy, dc, dr, dz))
  }) 
}

#Function for Monod estimation, quasi-steady state assumption
Monod_sts <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    #dx = mu_max_sts*x*(R_0 + 1/rho_sts*R_0 - 1/rho_sts*x)/(R_0 + 1/rho_sts*R_0 - 1/rho_sts*x + 2.5e-05) #2.5e-05 is the Ks value
    dx = kappa_2*x*r/(r + (kappa_2+kappa_3)/kappa_1) #
    dy = 0
    dc = 0
    dr = -(kappa_2 + kappa_3)*x*r/(r + (kappa_2+kappa_3)/kappa_1) #
    dz = 0
    list(c(dx, dy, dc, dr, dz))
  }) 
}

##===================================
## Fitted with logistic model #
##===================================
## numeric solution 
## ODEs system
parms_init <- c(mu_max = 0.5, rho = 0.5)
times <- TT

## model cost,
ModelCost2 <- function(P) {
  out <- ode(y = state, func = logist, parms = P, times = TT)
  model = out
  return(modCost(out, Data_Cells)) # object of class modCost
}

Fit <- modFit(f = ModelCost2, p = parms_init, lower = rep(0, 2),
              upper = c(5, 1))

out <- ode(y = state, func = logist, parms = Fit$par,
           times = Time_Eval)

plot(Data_Cells, xlim = c(0, 400), pch = 21, bg = alpha("green", 0.4), col = alpha("green", 0.4))
lines(out, col = "red", lty = 2)
summary(Fit)

##===================================
## Fitted with Monod model #
##===================================
## numeric solution 
## ODEs system

param_init = Fit$par
kappa_2 = param_init[1]
kappa_3 = kappa_2/0.5 - kappa_2[1]
kappa_1 = 1e3
parameters_init <- c(kappa_1 = as.numeric(kappa_1), kappa_2 = as.numeric(kappa_2), kappa_3 = as.numeric(kappa_3))
times <- TT

## model cost,
ModelCost3 <- function(P){
  out <- ode(y = state_1, func = Monod, parms = P, times = TT, atol = 1e-11, rtol = 1e-10)
  model = out
  model[,2] = model[,6]
  model = model[,1:2]
  return(modCost(model, Data_Cells)) # object of class modCost
}

Fit2 <- modFit(f = ModelCost3, p = parameters_init, lower = rep(0, 3), upper = c(1e10,50,50))

out2 <- ode(y = state_1, func = Monod, parms = Fit2$par,
           times = Time_Eval, atol = 1e-11, rtol = 1e-10)

lines(Time_Eval, out2[,2], col = "blue", lty = 2)

summary(Fit2)

##===================================
## Fitted with Monod model and steady state assumption#
##===================================
## numeric solution 
## ODEs system

param_init = Fit$par
#mu_max_sts = param_init[1]
#rho_sts = param_init[2]
kappa_2 = param_init[1]
kappa_3 = 2.666667e-05 #kappa_2/0.5 - kappa_2[1]
kappa_1 = 7.273933e-15 #1e2
parameters_init <- c(kappa_1 = as.numeric(kappa_1), kappa_2 = as.numeric(kappa_2), kappa_3 = as.numeric(kappa_3))
#parameters_init <- c(mu_max_sts = as.numeric(mu_max_sts), rho_sts = as.numeric(rho_sts))
times <- TT

## model cost,
ModelCost_sts <- function(P) {
  out <- ode(y = state_sts, func = Monod_sts, parms = P, times = TT, atol = 1e-11, rtol = 1e-10)
  model = out
  model = model[,1:2]
  return(modCost(model, Data_Cells)) # object of class modCost
}

Fit_sts <- modFit(f = ModelCost_sts, p = parameters_init, lower = rep(0, 3), upper = c(1e-10,50,150))

out_sts <- ode(y = state_sts, func = Monod_sts, parms = Fit_sts$par,
            times = Time_Eval, atol = 1e-11, rtol = 1e-10)

lines(Time_Eval, out_sts[,2], col = "purple", lty = 2)

summary(Fit_sts)

