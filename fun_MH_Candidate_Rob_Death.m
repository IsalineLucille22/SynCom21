function [accept_temp, W_n_thresh, W_n_Cov_Mat, Loop_Mat_Temp_cand, Loop_Mat_Temp, Threshold_Loop_cand, Threshold_Loop, Struct_Loop_Mat, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, p, n_loop, theta_T_Loop, theta_LT_Loop, nu, Cov_Mat, Cov_Mat_LT, LT_Loop_cand, LT_Loop] = fun_MH_Candidate_Rob_Death(accept_temp, k, W_n_thresh, W_n_Cov_Mat, Loop_Mat_Temp_cand, Loop_Mat_Temp, Threshold_Loop_cand, Threshold_Loop, Struct_Loop_Mat, p, n_loop, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, S_consumer, S_consumed, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, theta_T_Loop, theta_LT_Loop, nu, Cov_Mat, Cov_Mat_LT, max_val, nb_rep, LT_Loop_cand, LT_Loop, Weight_species, weight_day_seven, Temp)
%Function Metropolis Hasting. Compute the energy of the new candidate and
%compare it to the present state.      
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Energy determination %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % X = 1e3*X; Measured_Abund = 1e3*Measured_Abund;
    std_weight = ones(1, length(X(1,:)));
    errors = Quantile_Transformer(X, Measured_Abund, std_weight);
    errors_rel = Quantile_Transformer(StackPlotTot, StackPlot_Meas, std_weight);
    mean_measured_abund = mean(Measured_Abund(:, end, :), 3);
    mean_stackplot_meas = mean(StackPlot_Meas(:, end, :), 3);
    var_obs_abs = sum(mean_measured_abund.^2) + 1e-10;  % Regularization term added
    var_obs_rel = sum(mean_stackplot_meas.^2) + 1e-6;  % Regularization term added
    error_abs_sum = sum((errors.^2)./var_obs_abs, 'all')/nb_rep;
    error_rel_sum = sum((errors_rel.^2)./var_obs_rel, 'all')/nb_rep;
    alpha = 0.5;  % Weight between absolute and relative error components
    energie = alpha*error_abs_sum + (1 - alpha)*error_rel_sum;
    scaling_factor = max(mean_measured_abund.^2) + 1e-16;  % Prevent division by zero
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

    ratio = min(exp(-1*beta*delta_nrj/Temp), 1);
    % disp(ratio)
    if rand_val < ratio 
        Loop_Mat_Temp = Loop_Mat_Temp_cand;
        Threshold_Loop = Threshold_Loop_cand; %Uncomment it
        LT_Loop = LT_Loop_cand;
        acceptance_ratio_tot(p) = 1;
        acceptance_ratio_loop(n_loop) = 1;
        accept_temp(k) = 1;
        % Loop_Mat_Temp_mat = [Loop_Mat_Temp_mat, reshape(Loop_Mat_Temp,[],1)]; %Uncomment only if nb resources is fixed
    else
        liste_nrj_tot(p) = liste_nrj_tot(p - 1);
    end
    ratio = sum(accept_temp)/k;
    
    % Struct_Loop_Mat{n_loop + 1} = Loop_Mat_Temp;
    % n_loop = n_loop + 1;
    % p = p + 1;

    % Modification matrix coefficients
    nu_2 = min(1, 2*nu);
    N = S_consumer;
    Loop_Mat_Temp = diag(Loop_Mat_Temp);
    Cov_Mat = Cov_Mat*(eye(N) + nu_2*(accept_temp(k) - 0.234)*(W_n_Cov_Mat*W_n_Cov_Mat')/norm(W_n_Cov_Mat)^2)*Cov_Mat';%Cov_Mat*(eye(N) + nu_2*(ratio - 0.234)*(W_n_Cov_Mat*W_n_Cov_Mat')/norm(W_n_Cov_Mat)^2)*Cov_Mat';
    Cov_Mat = (Cov_Mat + Cov_Mat')/2;
    Cov_Mat = chol(Cov_Mat + 1e-6 * eye(N), 'lower');
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
        % if all(W_n_Cov_Mat == 0)
        %     c = 10;
        % end
    end

    % Modification thresholds
    theta_T_Loop = theta_T_Loop*(eye(S_consumer) + nu_2*(ratio - 0.234)*( W_n_thresh* W_n_thresh')/norm( W_n_thresh)^2)*theta_T_Loop';
    theta_T_Loop = chol(theta_T_Loop,'lower');
    W_n_thresh = normrnd(zeros(N,1), ones(N,1));
    prod_W_n = theta_T_Loop* W_n_thresh;
    temp = Threshold_Loop + prod_W_n;
    Threshold_Loop_cand = max(temp, 0); 
    W_n_thresh(temp<0) = prod_W_n(temp<0) - temp(temp<0);
    
    % Modification lag times
    LT_Loop_cand = LT_Loop;
    rand_inidice = rand(S_consumer, S_consumed) > 0.98;
    LT_Loop_cand(rand_inidice) = max(LT_Loop(rand_inidice) + normrnd(0, theta_LT_Loop, sum(sum(rand_inidice)), 1),0);
    logtheta = log(theta_LT_Loop) + nu*(ratio - 0.234);%
    theta_LT_Loop = exp(logtheta);%theta_LT_Loop(rand_inidice) = exp(logtheta);
end