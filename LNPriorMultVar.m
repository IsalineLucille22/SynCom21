function [delta, delta_test] = LNPriorMultVar(x_cand, x_obs, mu_Likelihood, var_Likelihood)
x_vect_cand = reshape(x_cand,1,[]);
non_zeros_val_1 = x_vect_cand > 0;
x_vect_obs = reshape(x_obs,1,[]);
non_zeros_val_2 = x_vect_obs > 0;
non_zeros_val = logical(non_zeros_val_1.*non_zeros_val_2);
x_vect_cand = x_vect_cand(non_zeros_val)';
x_vect_obs = x_vect_obs(non_zeros_val)';
mean_param_LN = mu_Likelihood(non_zeros_val)';
% n = length(x_vect);
% var_param_LN = diag(var_Likelihood(non_zeros_val));
% L = (2*pi)^(-n/2)*det(var_param_LN)^(-1/2)*prod(1./x_vect)*exp(-1/2*(log(x_vect) - mean_param_LN)'*(var_param_LN\(log(x_vect) - mean_param_LN)));
var_param_LN = var_Likelihood(non_zeros_val)';
L_cand = prod(lognpdf(x_vect_cand, mean_param_LN, var_param_LN));%prod(unifpdf(x_vect_cand, mean_param_LN, var_param_LN));%
L_obs = prod(lognpdf(x_vect_obs, mean_param_LN, var_param_LN));%prod(unifpdf(x_vect_obs, mean_param_LN, var_param_LN));%
delta = L_cand/L_obs;
delta_test = delta;
if L_obs == 0 && L_cand == 0
    delta = 1;
    delta_test = delta;
end
if isnan(delta)
    c = 10;
end
end