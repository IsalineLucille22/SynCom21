library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
# install.packages("Hmisc")
# library(Hmisc)
library(FME)
library(deSolve)

Dist_Euclid_mat <- function(mat){
  n = length(mat[,1]) 
  m = length(mat[1,]) 
  vect_sum_col = colSums(mat)
  mat_fin = diag(rep(0, n))
  for(j in 1:n){
    vect_res = rep(NA, j - 1)
    vect_res = rep(0,n)
    for(l in j:n){
      res = sum((mat[j,]/sum(mat[j,]) - mat[l,]/sum(mat[l,]))^2)
      res = sqrt(sum(res))
      res = res + (1 - cor(mat[j,], mat[l,])) + as.numeric((1 - mat[j,]%*%mat[l,]/(norm(as.numeric(mat[j,]), "2")*norm(as.numeric(mat[l,]), "2"))))
      vect_res[l] = res
    }
    mat_fin[j, ] = vect_res
  }
  mat_fin[lower.tri(mat_fin)] = t(mat_fin)[lower.tri(mat_fin)]
  diag(mat_fin) = 0
  return(mat_fin)
}

#Main function without executing the simulated annealing
K_Center_Clust<-function(Dist_Mat, T_max, nb_centers, centers){#Seeds of each cluster into the observations
  nb_obs = dim(Dist_Mat)[1]
  centers = centers[order(centers)]
  L = rep(0, nb_obs)
  wss = rep(0, nb_centers)
  ind_glob_mean = which(rowSums(Dist_Mat) == min(rowSums(Dist_Mat)))
  if(length(ind_glob_mean) > 1){
    ind_glob_mean = sample(ind_glob_mean,1)
  }
  for(i in 1:T_max){
    bss = 0
    diff_nb_centers = nb_centers - length(centers) 
    if(diff_nb_centers > 0){ #Randomly add centers if the number of centers doesn't correspond to nb_centers 
      vect_centers = 1:nb_obs
      vect_centers = vect_centers[!vect_centers %in% centers] #Add a center which is not already in the center set
      add_centers = sample(vect_centers, diff_nb_centers)
      centers = c(centers, add_centers)
    }
    centers = centers[order(centers)]
    for(j in 1:nb_obs){
      vect_temp = Dist_Mat[j, ]
      min = min(vect_temp[centers])
      vect_temp[-centers] = Inf
      ind_temp = which(vect_temp[]==min, arr.ind = T)
      if(length(ind_temp)==1){ #Si est plus proche d'un seul centre, alors va dans le cluster de ce centre
        L[j] = ind_temp
      }
      else{
        if(length(which(ind_temp == j)) == 1){ # Si une observation est un centre, alors va dans le cluster dont il est le centre
          L[j] = j
        }
        else{
          L[j]=sample(ind_temp,1) #Sinon choix aléatoire d'un cluster 
        }
      }
    }
    centers = unique(L)
    L_temp = L
    centers_temp = centers
    for(k in 1:nb_centers){
      ind = which(L==centers[k])
      mat_temp = as.matrix(Dist_Mat[ind,ind])
      vect_sum_dist = rowSums(mat_temp)
      ind_min = which(vect_sum_dist == min(vect_sum_dist)) #Calcul du nouveau centre du cluster
      if(length(ind_min)>1){ #If there are two centers. Two seeds with the same minimum distance, then randomly choose one of these two seeds.
        ind_min = sample(ind_min, 1)
      }
      ind_min = ind[ind_min]
      sum_min = min(vect_sum_dist)
      wss[k] = sum(Dist_Mat[ind_min,ind]^2) #Sum within the cluster
      centers_temp[k] = ind_min
      L_temp[ind] = ind_min
      bss = bss + length(ind)*Dist_Mat[ind_min,ind_glob_mean]^2
    }
    L = L_temp
    centers = centers_temp
    centers = unique(centers)
  }
  bss = 1/(nb_centers - 1)*bss
  wss = (1/(nb_obs - nb_centers))*sum(wss)
  wss = sum(wss)
  if(wss == 0){
    Var_var = 10000
  }
  else{
    Var_var = bss/wss
  }
  return(c(L,wss)) #What is displayed 
}

#Recuit_Simule_n involves the function K_Center_Clust
Recuit_Simule_n <- function(n_clust, Dist_Mat, n_iter, k_max){
  T_0 = 10#1000
  nb_obs = dim(Dist_Mat)[1]
  centers_1 = sample(1:nb_obs, n_clust)
  res = rep(0, k_max)
  Threshold = rep(0, k_max)
  Temp = T_0;
  val_comp = res[1]
  alpha = 0
  n_not_change = length(centers_1) - 2
  diff = n_clust - n_not_change
  for(i in 1:k_max){
    vect_centers = c(1:nb_obs)
    res_1 = K_Center_Clust(Dist_Mat,n_iter,n_clust, centers_1)
    L_1 = res_1[1:nb_obs]
    centers_1 = unique(L_1)
    E_1 = res_1[(nb_obs+1)]
    sample_kept = sample(centers_1, n_not_change)
    vect_centers = vect_centers[!vect_centers %in% sample_kept]
    add_centers = sample(as.numeric(vect_centers), diff)
    centers_2 = as.numeric(c(sample_kept, add_centers))
    res_2 = K_Center_Clust(Dist_Mat,n_iter,n_clust, centers_2)
    L_2 = res_2[1:nb_obs]
    E_2 = res_2[(nb_obs+1)]
    delta_E = -(E_2 - E_1)
    Threshold[i] = min(exp(delta_E/Temp),1)
    if(runif(1) <= exp(delta_E/Temp)){#Look at this more in detail. Put the log (exp(delta_E/Temp) >= runif(1))
      E_1 = E_2
      L_1 = L_2
      centers_1 = unique(L_2)
      alpha = alpha + 1 #We keep the candidate
    }
    if(E_1 != res[i]){
      val_comp = res[i]
    }
    res[i + 1] = E_1
    Temp = 0.9999*Temp #T_0/log(i + 1)
  }
  print(c(abs(res[i]-val_comp), alpha/k_max))
  plot(1:k_max, Threshold)
  return(c(L_1, E_1))
  #return(res[k_max])
}

Recuit_Simule <- function(vect_n_clust, D_Mat, n_iter, k_max){ #n_clust doit être >=2
  m = length(vect_n_clust)
  mat_res = rep(0, m)
  for(i in 1:m){
    res = Recuit_Simule_n(vect_n_clust[i], D_Mat, n_iter, k_max)
    mat_res[i] = res
  }
  return(mat_res)
}

#Silhouette method when we don't know the exact number of cluster
Silhouette_meth_n <- function(D_Mat, T_max, nb_centers){
  m = nb_centers
  n = length(D_Mat[1,]) # Nombre total d'observations 
  centers = sample(1:n, nb_centers)
  L = Recuit_Simule_n(nb_centers, D_Mat, 1, T_max)
  L = L[1:n]
  Clust_renorm = L
  num_clust = unique(Clust_renorm)
  s_sil = rep(0, n) # Vecteur silhouette de chaque élément
  S_sil = 0
  for(i in 1:n){
    ind_mc = which(Clust_renorm == Clust_renorm[i]) #Find the observations from the same cluster
    clust_diff = num_clust[which(num_clust != Clust_renorm[i])] # Indique les différents clusters
    if(length(ind_mc) == 1){
      a = 0
    }
    else{
      a = sum(D_Mat[i, ind_mc]) #Somme des distances aux éléments d'un même cluster
      a = (1/(length(ind_mc) - 1))*a
    }
    b = Inf
    for(j in 1:(m-1)){ #length(clust_diff)) = m-1
      ind_temp = which(Clust_renorm == clust_diff[j])
      nb_elem_in_clust_temp = length(ind_temp)
      b_temp = sum(D_Mat[i, ind_temp])
      b_temp = (1/nb_elem_in_clust_temp)*b_temp
      if(b_temp < b){
        b = b_temp
      }
    }
    if(a == 0){
      res = 0
      #res = (b-a)/(max(a,b))
    }
    else{
      res = (b-a)/(max(a,b))
    }
    s_sil[i] = res
  }
  for(i in 1:m){
    ind_temp =  which(Clust_renorm == num_clust[i])
    nb_elem = length(ind_temp)
    res = sum(s_sil[ind_temp])
    S_sil = S_sil + res/nb_elem
  }
  S_sil = (1/m)*S_sil
  return(S_sil)
}

Silhouette_meth <- function(D_Mat, T_max, vect_centers){
  n = length(vect_centers)
  S_sil = rep(0,n)
  for(i in 1:n){
    S_sil[i] = Silhouette_meth_n(D_Mat, T_max, vect_centers[i])
  }
  return(S_sil)
}


SSE<-function(vect_n_clust, Mat_obs, D_Mat, n_iter, k_max){
  m = length(vect_n_clust)
  mat_res = rep(0, m)
  nb_obs = dim(Dist_Mat)[1]
  for(i in 1:m){
    res = Recuit_Simule_n(vect_n_clust[i], D_Mat, n_iter, k_max)
    L_temps = res[1:nb_obs]
    centers_temps = unique(L_temps) #Center among observations, not the euclidian center
    sse_nb_cluster = 0 #Sum Standard Error for a given number of clusters
    for(j in 1:vect_n_clust[i]){
      sse = 0 #Sum Standard Error for a given cluster
      ind_temp = which(L_temps == centers_temps[j])
      if(length(ind_temp) > 1){
        true_center = colSums(Mat_obs[ind_temp, ])/length(ind_temp)
        sse = (Mat_obs[ind_temp, ] - matrix(rep(true_center, length(ind_temp)), nrow = length(ind_temp), byrow = TRUE))^2
        sse = sum(rowSums(sse))
      } 
      else{
        true_center = Mat_obs[ind_temp, ]
        sse = 0
      }
      sse_nb_cluster = sse_nb_cluster + sse
    }
    mat_res[i] = sse_nb_cluster
  }
  return(mat_res)
}

Renorm_Clust <- function(vect){
  Ind_Clust = unique(vect)
  for(i in 1:length(Ind_Clust)){
    ind = which(vect == Ind_Clust[i])
    vect[ind] = 1000 + i
  }
  vect = vect - rep(1000, length(vect))
  return(vect)
}


setwd("/Users/iguex/Documents/CoCulture_Soil/Data")
Data = read_excel("MergedData.xlsx", sheet = 11, range = "A1:CU21")
Parameters_set = Data;
Mat_list = data.matrix(Parameters_set[,2:98])
# Dist_Mat = as.matrix(dist(Mat_list))
Dist_Mat = Dist_Euclid_mat(Mat_list)

nb_clust = c(3:19)
res_Sil = Silhouette_meth(Dist_Mat, 50000, nb_clust)#270000
plot(nb_clust, res_Sil, type = "l")
abline(v = 13, col = "gray", lty = 2)
abline(v = 12, col = "gray", lty = 2)

test_new_Euclid_R_12_2 = Recuit_Simule_n(6, Dist_Mat, 1, 25000)#300000
res_12_2 = Renorm_Clust(test_new_Euclid_R_12_2[1:20])#211
print(res_12_2)
print(test_new_Euclid_R_12_2[21])

