function [accept_temp, W_n_Cov_Mat, Loop_Mat_Temp_cand, Loop_Mat_Temp, Struct_Loop_Mat, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, p, n_loop, nu, Cov_Mat, param_chain] = fun_MH_Candidate_Rob_Death(accept_temp, k, W_n_Cov_Mat, Loop_Mat_Temp_cand, Loop_Mat_Temp, Struct_Loop_Mat, p, n_loop, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, S_consumer, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, nu, Cov_Mat, max_val, nb_rep, Weight_species, weight_day_seven, Temp, param_chain, E_grad)
%Function Metropolis Hasting. Compute the energy of the new candidate and
%compare it to the present state.      
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Energy determination %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    std_weight = ones(1, length(X(1,:)));
    errors = Quantile_Transformer(X, Measured_Abund, std_weight);
    mean_measured_abund = mean(Measured_Abund(:, end, :), 3);
    var_obs_abs = sum(mean_measured_abund.^2);  % Regularization term added
    error_abs_sum = sum((errors.^2)./var_obs_abs, 'all')/nb_rep;
    energie = error_abs_sum;
    scaling_factor = 1e-03;%max(mean_measured_abund.^2) + 1e-16;  % Prevent division by zero
    energie = energie/scaling_factor;
    lambda_L2 = 1e-03;
    energie = energie + lambda_L2*sum(Loop_Mat_Temp_cand(:).^2)/numel(Loop_Mat_Temp_cand);
    energie = max(energie, 1e-10);

    rand_val = unifrnd(0,1);
    liste_nrj_loop(n_loop) = energie; 
    liste_nrj_tot(p) = energie;

    if p > 1
        delta_nrj = liste_nrj_tot(p) - liste_nrj_tot(p-1); 
    else
        delta_nrj = 0;
    end

    % prior_CF = exp(-lambda_CF * sum((Loop_Mat_Temp_cand(:) - mu_CF).^2));
    % ratio = min(exp(-beta*delta_nrj/Temp)*prior_CF, 1);
    ratio = min(exp(-1*beta*delta_nrj/Temp), 1);
    disp(ratio)
    if rand_val < ratio 
        Loop_Mat_Temp = Loop_Mat_Temp_cand;
        acceptance_ratio_tot(p) = 1;
        acceptance_ratio_loop(n_loop) = 1;
        accept_temp(k) = 1;
    else
        liste_nrj_tot(p) = liste_nrj_tot(p - 1);
    end
    param_chain(:, :, p) = Loop_Mat_Temp_cand;
    
    % Modification matrix coefficients
    nu_2 = min(1, 2*nu);
    N = S_consumer;
    Loop_Mat_Temp = diag(Loop_Mat_Temp);
    Cov_Mat = Cov_Mat*(eye(N) + nu_2*(accept_temp(k) - 0.44)*(W_n_Cov_Mat*W_n_Cov_Mat')/norm(W_n_Cov_Mat)^2)*Cov_Mat';
    Cov_Mat = (Cov_Mat + Cov_Mat')/2;
    Cov_Mat = chol(Cov_Mat + 1e-4*eye(N), 'lower');
    W_n_Cov_Mat = zeros(N,1);
    while all(W_n_Cov_Mat == 0)
        W_n_Cov_Mat = normrnd(zeros(N,1), 0.1*ones(N,1));
        prod_W_n = Cov_Mat*W_n_Cov_Mat;
        temp = Loop_Mat_Temp + prod_W_n;
        Loop_Mat_Temp_cand = max(temp, 0); %Cov_Mat*W_n only change a fraction of indices
        W_n_Cov_Mat(temp<0) = prod_W_n(temp<0) - temp(temp<0);
        % W_n_Cov_Mat(temp < 0) = 0.5*(prod(temp < 0) - temp(temp < 0));
        High_val = find(Loop_Mat_Temp_cand > max_val);
        Loop_Mat_Temp_cand(High_val) = Loop_Mat_Temp(High_val); %Upper boundary
        Loop_Mat_Temp = diag(Loop_Mat_Temp);
        Loop_Mat_Temp_cand = diag(Loop_Mat_Temp_cand);
        W_n_Cov_Mat(High_val) = 0;
    end
end