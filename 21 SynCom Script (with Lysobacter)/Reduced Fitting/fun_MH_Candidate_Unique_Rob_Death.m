function [accept_temp, W_n_Cov_Mat, Loop_Mat_Temp_cand, Loop_Mat_Temp, Struct_Loop_Mat, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, p, n_loop, nu, Cov_Mat, param_chain] = fun_MH_Candidate_Unique_Rob_Death(accept_temp, k, W_n_Cov_Mat, Loop_Mat_Temp_cand, Loop_Mat_Temp, Struct_Loop_Mat, p, n_loop, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, S_consumer, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, nu, Cov_Mat, max_val, nb_rep, Weight_species, weight_day_seven, Temp, param_chain, E_grad)
%Function Metropolis Hasting. Compute the energy of the new candidate and
%compare it to the present state.   
% =========================================================
% Metropolis–Hastings for ONE parameter
% =========================================================
% 
% % Initialize
% theta = theta_init;        % current parameter value
% sigma = 0.1;               % proposal std (step size)
% alpha_target = 0.44;       % target acceptance for 1D
% nu = 0.01;                 % adaptation step size
% beta = 1;
% Temp = 1;
% lambda_L2 = 1e-3;

% sigma = 5e-1;               % proposal std (step size)
% nu_in = 1;%0.01;
% alpha_target = 0.44;
errors = Quantile_Transformer(X, Measured_Abund, ones(1, length(X(1,:))));
mean_measured_abund = mean(Measured_Abund(:, end, :), 3);
var_obs_abs = sum(mean_measured_abund.^2);
error_abs_sum = sum((errors.^2)./var_obs_abs, 'all')/nb_rep;
scaling_factor = 1e-3;%max(mean_measured_abund.^2) + 1e-16;  % Prevent division by zero
error_abs_sum = error_abs_sum/scaling_factor;
lambda_L2 = 1e-03;
energie = error_abs_sum + lambda_L2*Loop_Mat_Temp_cand^2;
energie = max(energie, 1e-10);

rand_val = unifrnd(0,1);
liste_nrj_loop(n_loop) = energie; 
liste_nrj_tot(p) = energie;

if p > 1
    delta_nrj = liste_nrj_tot(p) - liste_nrj_tot(p-1); 
else
    delta_nrj = 0;
end

log_alpha = -1e01*beta*delta_nrj/Temp;
ratio = min(1, exp(min(0, log_alpha)));
% disp(ratio)
if rand_val < ratio 
    Loop_Mat_Temp = Loop_Mat_Temp_cand;
    acceptance_ratio_tot(p) = 1;
    acceptance_ratio_loop(n_loop) = 1;
    accept_temp(k) = 1;
else
    liste_nrj_tot(p) = liste_nrj_tot(p - 1);
end
param_chain(:, :, p) = Loop_Mat_Temp_cand;

% --- Adaptive update of step size (Robbins–Monro style) ---
% sigma = sigma*exp(nu_in*(accept_temp(k) - alpha_target));
% added = normrnd(0, 0.1);
% Loop_Mat_Temp_cand = max(Loop_Mat_Temp + sigma*added, 0); %max(Loop_Mat_Temp + sigma*randn(), 0); 


nu_2 = min(1, 2*nu);
Cov_Mat = Cov_Mat * sqrt(1 + nu_2*(accept_temp(k) - 0.44));
W_n_Cov_Mat = randn();              % scalar standard normal
prod_W_n = Cov_Mat * W_n_Cov_Mat;   % scalar step
temp = Loop_Mat_Temp + prod_W_n;    % proposed new value
Loop_Mat_Temp_cand = max(temp, 1e-06);
if Loop_Mat_Temp_cand > max_val
    Loop_Mat_Temp_cand = Loop_Mat_Temp;  % reject if outside range
end
end