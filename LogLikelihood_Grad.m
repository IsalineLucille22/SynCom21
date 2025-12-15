function grad_vect = LogLikelihood_Grad(errors, errors_2, var_obs_abs, var_obs_rel)
grad_vect = 4*sum(errors)/(var_obs_abs * nb_rep*1e06) + 2*2/nb_rep*sum(errors_2)/var_obs_rel; %Look at sum
end