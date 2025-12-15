clear;
close all;

%Determine the resources consumption and death matrix on the monocultures
%Save or Not
save_data = 1; %1 if save, 0 otherwise %Do not save for now
% Create a structure or cell array to hold the data
data_to_save_Mono = struct(); %To avoid overwriting, change the name of the save variable at the last line of the script. Name it accoridng to this structure name.
addpath('/Users/pret_helpdesk/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
addpath('/Users/pret_helpdesk/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')
addpath('/Users/pret_helpdesk/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/21 SynCom Script (with Lysobacter)/Data')


%Parameters drawn from SnyCom20 fitting
Name_file =  'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Newv2';
load(strcat('Data/', Name_file));
% parfor nb_iter_tot = 1:20   
for nb_iter_tot = 1:1

    %Loadind data
    Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
    Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');
    Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','142:163', 'Format','auto');
    Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 11, 'Range','1:22', 'Format','auto');%Monoculture growth in soil
    Data_Evol_test = Data_Evol;
    Time_step = [0 12 72 168 321]; %Measured time step in hours 
    tspan = [0, max(Time_step)]; %Time interval in hours
    S = height(Data_Evol);
    %Load previously and assumed correct lag times and thresholds
    load(strcat('Data/', Name_file));
    data_to_save = data_to_save_SynCom21_Soil_only_Lyso_Parfor_New;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_New;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_7;%data_to_save_Inter_v2;%data_no_CF;
    % Fitted parameters
    ind_sim = 3; %We can choose any ind_sim because threshold values are constant
    Threshold_CF = data_to_save(ind_sim).Threshold_CF;
    Threshold_death = data_to_save(ind_sim).Threshold_death; 
    Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons;
    Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; %0.001*eye(S,S);%
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix;
    Pred_Mat_Lyso = zeros(S, S);
    yield_Pred = 0.2;

    %Turn 7 into 5
    yield_vect = table2array(Parameters_set(1:S,7)); %Change it according to the desired mu_max
    name = string(table2array(Parameters_set(1:S,1)));%string(table2array(table(2:22,3)));
    
    mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8)); %Log-normal Parameters_Senka_mu_max(:,7:8)
    mu_max_vect = mu_max_dist(:,1);
    Res_Percentage = 1;%table2array(Parameters_Senka_mu_max(:,9)); %1 Assuming each resource can be consumed %Percentage of resources that can be consumed
    Parameters_Senka_Lag_time = table2array(Parameters_Senka_Lag_time(:,7:8)); %table2array(Parameters_Senka_Lag_time(:,2:3));
    mean_param_LN = mu_max_dist(:,1);
    var_param_LN = mu_max_dist(:,2);
    
    Data_Evol_temp = table2array(Data_Evol(:, 2:end));
    nb_obs = length(Data_Evol_temp(1,:));
    nb_time_step = length(Time_step);
    nb_rep = nb_obs/nb_time_step;
    Measured_Abund = zeros(length(mu_max_vect), nb_time_step, nb_rep); %Number species, number times, number replicates.
    for i = 1:nb_rep
        Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
    end
    
    mean_y_0 = mean(Measured_Abund(:,1,:), 3);%table2array(Data_Evol(1:20, 2));
    std_y_0 = std(Measured_Abund(:,1,:), 1, 3);  
    StackPlot_Meas = Measured_Abund./sum(Measured_Abund);

    %%%%%=== PRECOMPUTE DERIVATIVE PSEUDO-DATA === 
    % === Build nonuniform pseudo-grid proportional to observed intervals ===
    t_obs = Time_step';
    G_target = nb_time_step*4;   % about 4× more pseudo-points than observations
    t_grid = t_obs(1);
    dt = diff(t_obs); %Difference between two time points
    Tsum = sum(dt);
    for k = 1:length(dt)
        % number of internal points in this interval proportional to its length
        gk = max(1, round(G_target*dt(k)/Tsum));
        tk = linspace(t_obs(k), t_obs(k+1), gk + 2);
        t_grid = [t_grid; tk(2:end-1)'];   % internal points only
    end
    t_grid = unique([t_grid; t_obs]);      % ensure observed times included


    %Spline reconstruction of pseudo-data. Around 3-5*nb_time_step pseudo
    %points
    S = size(Measured_Abund, 1);
    MeanTraj = mean(Measured_Abund, 3);  % S × T matrix
    
    dx_hat = zeros(S, numel(t_grid));
    for i = 1:S
        yi = squeeze(MeanTraj(i, :))';
        fitobj = fit(t_obs, yi, 'smoothingspline');
        dx_hat(i, :) = differentiate(fitobj, t_grid);
    end
    
    % Estimate derivative variance per species
    var_dx = var(diff(Measured_Abund,1,2), 0, [2 3]) + 1e-8;
    
    % Save for reuse in the MH function
    save('pseudo_grad_data.mat', 'dx_hat', 't_grid', 'var_dx');
    %%%%%


    liste_nrj_tot_temp = 1e15;
    %Initialization of the model parameters fixed for all replicates 
    t_0 = 0; %Time 0
    nb_Res = 12; %Number of resources (group of resources)
    CrossFeed_Mat = lognrnd((log(0.8) + mu_max_dist(:,1)).*ones(S,S), 0*mu_max_dist(:,2).*ones(S,S));%.val val% of zeros. ConsumerxProducer
    death_rate = -[0;-0.000917642333192334;0;0;-0.000631406110652123; -0.000354932;-0.00122911113259887;...
        -0.00142098782776364;-0.00168390810416051;-0.000942517121033203;0;-0.000720301703810259;0;...
        0;-0.000192751969464189;-0.00254373051637326;-0.00163181578562513;0;-0.00132840048140091;-0.00185805721783005;-0.00623013508583033];
    death_rate_2 = [1e-05; 6.395e-04; 3.5730e-04; 7.1931e-04;  8.2640e-05; 3.54932e-04; 6.261921e-04;  4.6248e-05; 2.6131e-04; 3.4113e-04;...
    9.4029e-04; 9.9890e-05; 3.8416e-04; 1e-05; 6.0388e-06; 5.8124e-04; 2.1584e-04; 6.5872e-04; 2.6565e-04; 1.01239e-03; 1.30413e-03];
    death_rate = mean([death_rate, death_rate_2], 2);
    CrossFeed_Mat_Temp = zeros(S,S); %No cross-feeding in Monocultures
    death_rate = max(death_rate, 7e-05);
    
    % Create a resource matrix
    %Load parameters
    R_0 = 0.5;
    % x = rand(S,nb_Res);
    % Res_Percentage_old = Res_Percentage;
    % Res_Percentage = min(max(Res_Percentage, 0.24), 0.24);
    % rand_indice = x > Res_Percentage;
    % Resource_Matrix = ones(S,nb_Res);
    % Resource_Matrix(rand_indice) = 0;
    % Resource_Matrix = Resource_Matrix.*lognrnd(-mu_max_dist(:,2).^2./2.*ones(S, nb_Res), mu_max_dist(:,2).*ones(S, nb_Res));
    R = R_0*15*10^(-3).*[0.6990 0.0581 0.0029 0.1477 4.0011e-05 0.0061 0.0233 0.0288 0.0008 0.0168 0.0070 0.0094]; %based on the observed values
    x = rand(S,nb_Res);
    rand_indice = x > 0.3; %Res_Percentage_Nb_Cons;
    Lag_time_Rand = normrnd(50, 10, S, nb_Res);%lognrnd(log(50), 1, S, nb_Res);
    Lag_time_Cons(rand_indice) = max(Lag_time_Cons(rand_indice), Lag_time_Rand(rand_indice));
    
    %Setting for Matlab ODE solver
    opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
    
    %Initialization of the colors
    
    colors = distinguishable_colors(S);
    
    num_fig = 1;
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0;
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];
    
    %Fit one unique threshold value
    Threshold_Pred = max(normrnd(1e-04, 0, S, 1),0);
    mat_y_0 = reshape(mat_y_0', 1, []);
    mat_y_0 = [mat_y_0 R];
    mu_max_vect_temp = exp(mu_max_vect + mu_max_dist(:,2).^2./2);
    yield = max(normrnd(yield_vect, 0.0), 0);
    kappa_mat = [2.5e+05*ones(S,1) mu_max_vect_temp (mu_max_vect_temp./yield - mu_max_vect_temp) 2.5e+05*ones(S,1)];
    kappa_mat_init = kappa_mat; 
    Mat_kappa_3_tot = repmat(kappa_mat(:,3), 1, S);%repmat((kappa_temp./yield' - kappa_temp), 1, S);%Kappa_3 rates defined by species not byproduct.
    Mat_kappa_3 = Mat_kappa_3_tot;
    Mat_kappa_3(CrossFeed_Mat_Temp == 0) = 0;Mat_kappa_3_init = Mat_kappa_3;
        
    beta = 0.5*1e2; 
    num_steps = 1000;
    nb_loops = 2; %Number of different loops (CF, death, resources, predation), 3 if predation loop is commented
         
    n = 1; r = 1; l = 1; p = 1; u = 1; fig_tot = 1; %Total energy
    max_tot_iter = 5;
    liste_nrj_CF = zeros(1, num_steps*max_tot_iter);
    liste_nrj_pred = zeros(1, num_steps*max_tot_iter);
    Struct_Death_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
    Struct_Death_Mat{1} = Death_Mat_Temp; Death_Mat_init = Death_Mat_Temp;
    Death_Mat_Temp_cand = Death_Mat_Temp; %Initilization candidat
    acceptance_ratio_death = zeros(1, num_steps*max_tot_iter);
    liste_nrj_res = zeros(1, num_steps*max_tot_iter);
    Struct_Res_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
    Struct_Res_Mat{1} = Resource_Matrix; Resource_Matrix_init = Resource_Matrix;
    Resource_Matrix_cand = Resource_Matrix; %Initilization candidat
    acceptance_ratio_res = zeros(1, num_steps*max_tot_iter);
    %%%
    liste_nrj_tot = zeros(1, nb_loops*num_steps*max_tot_iter);
    Lag_time_Cons_cand = Lag_time_Cons; Lag_time_Cons_init = Lag_time_Cons;
    mat_y_0_new = mat_y_0; 

    var_theta = 0.1;
    acceptance_ratio_temp = 10000;
    ratio_increase_loop = 1;
    N = nb_loops*num_steps*max_tot_iter;
    max_val = [10./mu_max_vect_temp 10*ones(S,1) 10./mu_max_vect_temp];
    Weight_species =  ones(S,8);
    weight_day_seven = 1;
    Death_Trace_Plot = []; 
    lag = 1;
    param_chain_2 = zeros(S, S, N);
    param_chain_1 = zeros(S, nb_Res, N);
    nb_Species = S;
    for iter_sp = 1:nb_Species
        S = 1;
        %Initial concentrations using a normal distribution
        mat_y_0 = mean_y_0(iter_sp);
        mat_y_0 = [mat_y_0 zeros(1,1) zeros(S,1)];
        
        %Fit one unique threshold value
        Threshold_Pred = max(normrnd(1e-04, 0, S, 1),0);
        mat_y_0 = reshape(mat_y_0', 1, []);
        mat_y_0 = [mat_y_0 R];
        % p = 1e10;
        while p < N
    
            %%%%%%%%%% Resource matrix loop %%%%%%%%
        
            Cov_Mat = eye(S*nb_Res);
            Cov_Mat_LT = 0.1*eye(S*nb_Res);
            W_n_LT_Res = normrnd(zeros(S*nb_Res,1), ones(S*nb_Res,1)); W_n_Cov_Mat_Res = normrnd(zeros(S*nb_Res,1), ones(S*nb_Res,1));
            accept_temp = zeros(1, num_steps);
            Temp = 500;
            for k = 1:(num_steps)
                for j = 1:(lag)
                    sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat(iter_sp, :), CrossFeed_Mat_Temp(iter_sp,:), Mat_kappa_3(iter_sp,:), Resource_Matrix_cand(iter_sp,:),...
                        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp(iter_sp,:), death_rate(iter_sp),...
                        Pred_Mat_Lyso(iter_sp,:), yield_Pred, S, Lag_time_Cons_cand(iter_sp,:), 0, nb_Res, 10, 10, mean_y_0(iter_sp)), tspan,  mat_y_0, opts_1);
                    z_temp = deval(sol, Time_step);
                    X = z_temp(mod(1:S*3,3) == 1, :); %Species'biomass
                    P = z_temp(mod(1:S*3,3) == 2, :); %Complexes'biomass
                    X = X + P;
                    nu = k^(-3/4); 
                    StackPlotTot = 1;
    
                    % === GRADIENT MATCHING PENALTY ===
                    lambda_grad = 1e-4;  % tune between 1e-2 and 1 
                    f_theta = zeros(S, numel(t_grid));
                    
                    for g = 1:numel(t_grid)
                        % Evaluate state from ODE solution
                        y_curr = deval(sol, t_grid(g));  % same as in your code below
                        % Compute derivative using your existing ODE RHS
                        dz_curr = fun_CF_Death_Lyso(t_grid(g), y_curr, kappa_mat(iter_sp, :), CrossFeed_Mat_Temp(iter_sp,:), Mat_kappa_3(iter_sp,:), Resource_Matrix_cand(iter_sp,:),...
                        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp(iter_sp,:), death_rate(iter_sp),...
                        Pred_Mat_Lyso(iter_sp,:), yield_Pred, S, Lag_time_Cons_cand(iter_sp,:), 0, nb_Res, 10, 10, mean_y_0(iter_sp));
                        % Extract only the dx part (species biomass derivative)
                        f_theta(:, g) = dz_curr(1:S);
                    end
                    % before computing diff_grad
                    dx_hat_i = dx_hat(iter_sp, :);        % 1 × G
                    f_theta_i = f_theta(1, :);            % 1 × G  (since S==1 inside iter_sp loop)
                    var_dx_i  = var_dx(iter_sp);          % scalar
                    
                    diff_grad = (dx_hat_i - f_theta_i).^2 ./ (var_dx_i + 1e-8);
                    % Normalize by grid length to keep scale comparable across settings
                    E_grad = lambda_grad*mean(diff_grad);
                    %%%
    
                    [accept_temp, W_n_LT_Res, W_n_Cov_Mat_Res, Resource_Matrix_cand(iter_sp,:), Resource_Matrix(iter_sp,:), Struct_Res_Mat, liste_nrj_tot, liste_nrj_res, ~, acceptance_ratio_res, p, l, nu, Cov_Mat, param_chain_1(iter_sp,:, :), Cov_Mat_LT, Lag_time_Cons_cand(iter_sp,:), Lag_time_Cons(iter_sp,:)] = ...
                        fun_MH_Candidate_Rob(accept_temp, k, W_n_LT_Res, W_n_Cov_Mat_Res, Resource_Matrix_cand(iter_sp,:), Resource_Matrix(iter_sp,:), Struct_Res_Mat, p, l, liste_nrj_tot, liste_nrj_res,...
                        0, acceptance_ratio_res, S, nb_Res, X, Measured_Abund(iter_sp,:, :), StackPlotTot, 1, beta, nu, Cov_Mat, Cov_Mat_LT, max_val(iter_sp, 3), nb_rep, Lag_time_Cons_cand(iter_sp,:), Lag_time_Cons(iter_sp,:), Weight_species, weight_day_seven, Temp, param_chain_1(iter_sp,:, :), E_grad);
                end
                Temp = Temp*0.99;
                % if mod(p, 50) == 0 && std(liste_nrj_tot(p - 49:p)) < 1e-8
                %     Temp = Temp;%Temp*0.5;
                % end
    
                Struct_Res_Mat{r + 1} = Resource_Matrix;
                l = l + 1;
                p = p + 1;
    
    
                if mod(p, num_steps) == 0
                    disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
                end
            end
            disp(sum(accept_temp)/(num_steps))
            burn_in = round(0.3*size(param_chain_1, 2));
            param_chain_post = param_chain_1(:, :, burn_in:fig_tot);
            
            % Compute mean, std, and CV for each parameter
            param_mean = mean(param_chain_post, 2);
            param_std  = std(param_chain_post, 0, 2);
            CV_3 = param_std./abs(param_mean);
            %%%
    
            %%%%%%%%%% Death loop %%%%%%%%
            
            Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp./mu_max_vect_temp;
            Cov_Mat = 0.5*eye(S);
            W_n_Cov_Mat_death = normrnd(zeros(S,1), ones(S,1)); 
            accept_temp = zeros(1, num_steps);
            Temp = 500;
            for k = 1:num_steps

                for j = 1:lag
                    sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat(iter_sp, :), CrossFeed_Mat_Temp(iter_sp,:), Mat_kappa_3(iter_sp,:), Resource_Matrix(iter_sp,:),...
                        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp_cand(iter_sp,:), death_rate(iter_sp),...
                        Pred_Mat_Lyso(iter_sp,:), yield_Pred, S, Lag_time_Cons(iter_sp,:), 0, nb_Res, 10, 10, mean_y_0(iter_sp)), tspan,  mat_y_0, opts_1);
                    z_temp = deval(sol, Time_step);
                    X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
                    P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
                    X = X + P;
                    nu = k^(-3/4); 
                    StackPlotTot = X./sum(X);

                    % === GRADIENT MATCHING PENALTY ===
                    lambda_grad = 1e-1;  % tune between 1e-2 and 1
                    f_theta = zeros(S, numel(t_grid));

                    for g = 1:numel(t_grid)
                        % Evaluate state from ODE solution
                        y_curr = deval(sol, t_grid(g));  % same as in your code below
                        % Compute derivative using your existing ODE RHS
                        dz_curr = fun_CF_Death_Lyso(t_grid(g), y_curr, kappa_mat(iter_sp, :), CrossFeed_Mat_Temp(iter_sp,:), Mat_kappa_3(iter_sp,:), Resource_Matrix(iter_sp,:),...
                        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp_cand(iter_sp,:), death_rate(iter_sp),...
                        Pred_Mat_Lyso(iter_sp,:), yield_Pred, S, Lag_time_Cons(iter_sp,:), 0, nb_Res, 10, 10, mean_y_0(iter_sp));
                        % Extract only the dx part (species biomass derivative)
                        f_theta(:, g) = dz_curr(1:S);
                    end

                    % before computing diff_grad
                    dx_hat_i = dx_hat(iter_sp, :);        % 1 × G
                    f_theta_i = f_theta(1, :);            % 1 × G  (since S==1 inside iter_sp loop)
                    var_dx_i  = var_dx(iter_sp);          % scalar

                    diff_grad = (dx_hat_i - f_theta_i).^2 ./ (var_dx_i + 1e-8);
                    % Normalize by grid length to keep scale comparable across settings
                    E_grad = lambda_grad*mean(diff_grad);
                    %%%

                    [accept_temp, W_n_Cov_Mat_death, Death_Mat_Temp_cand(iter_sp, iter_sp), Death_Mat_Temp(iter_sp, iter_sp), Struct_Death_Mat, liste_nrj_tot, liste_nrj_pred, ~, acceptance_ratio_death, p, r, nu, Cov_Mat, param_chain_2(iter_sp,:, :)] = ...
                        fun_MH_Candidate_Unique_Rob_Death(accept_temp, k, W_n_Cov_Mat_death, Death_Mat_Temp_cand(iter_sp, iter_sp), Death_Mat_Temp(iter_sp, iter_sp), Struct_Death_Mat, p, r, liste_nrj_tot, liste_nrj_pred,...
                        0, acceptance_ratio_death, S, X, Measured_Abund(iter_sp,:, :), StackPlotTot, 1, beta, nu, Cov_Mat, max_val(iter_sp, 2), nb_rep, Weight_species, weight_day_seven, Temp, param_chain_2(iter_sp,:, :), E_grad);

                    Death_Trace_Plot(:,:, r) = Death_Mat_Temp;
                    ind = eye(S,S);
                    Death_Mat_Temp_cand(~ind) = 0; %Changement pour diagonale zeros
                end
                % if mod(p, 50) == 0 && std(liste_nrj_tot(p - 49:p)) < 1e-8
                %     Temp = Temp*0.5;
                % end
                Temp = Temp*0.99;

                Struct_Death_Mat{r + 1} = Death_Mat_Temp;
                r = r + 1;
                p = p + 1;

                if mod(p, num_steps) == 0
                    disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
                end
            end
            disp(sum(accept_temp)/num_steps)
            % burn_in = round(0.3*size(param_chain_2, 2));
            % param_chain_post = param_chain_2(:, :, burn_in:r);

            % Compute mean, std, and CV for each parameter
            param_mean = mean(param_chain_post, 2);
            param_std  = std(param_chain_post, 0, 2);
            CV_2 = param_std./abs(param_mean);
            %%%%%%%%%%%%%%%%
        
            Resource_Matrix_cand = Resource_Matrix; %Re-initialization if the last candidate is not accepted. Otherwise issues with the resource loop.
            Lag_time_Cons_cand = Lag_time_Cons;
            liste_nrj_tot_temp = liste_nrj_tot(p - 1);
            num_steps = round(ratio_increase_loop*num_steps); 
            weight_day_seven = ratio_increase_loop*weight_day_seven;
            ratio_obs_sim = 1 - abs(X - mean(Measured_Abund,3))./(abs(X) + abs(mean(Measured_Abund,3))); %Similaritiy intex used as weights for variances
            Weight_species = mean(ratio_obs_sim, 2);
        end 
        
        %Comparison to observed data
        sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat(iter_sp, :), CrossFeed_Mat_Temp(iter_sp,:), Mat_kappa_3(iter_sp,:), Resource_Matrix(iter_sp,:),...
            Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp(iter_sp,:), death_rate(iter_sp),...
            Pred_Mat_Lyso(iter_sp,:), yield_Pred, S, Lag_time_Cons(iter_sp,:), 0, nb_Res, 10, 10, mean_y_0(iter_sp)), tspan,  mat_y_0, opts_1);
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        X = X + P;
                
        %%% Figure generation %%%
     
        h_1 = plot(Time_step, X, '-o', 'Color', colors(iter_sp,:), 'LineWidth', 1.5);%, Time_step, R_temp, 'o');
        axis([0 500 0 2.5e-03])
        % plot(Time_step, mean(Measured_Abund(iter_sp,:,:), 3), '--*', 'Color', colors(iter_sp,:), 'LineWidth', 1.5);%, Time_step, R_temp, 'o');
        % axis([0 500 0 2.5e-03])
        hold on 

        % num_fig = num_fig + 1; 
        p = 1;
    end
    for k = 1:nb_Species
        plot(Time_step, mean(Measured_Abund(k,:,:), 3), '--*', 'Color', colors(k,:), 'LineWidth', 1.5);
        hold on
    end
    axis([0 500 0 0.5e-03])
    legend(name, 'Orientation', 'vertical', 'Location', 'southeast')

    data_to_save_Mono(nb_iter_tot).Death_Mat_Temp = Death_Mat_Temp;
    data_to_save_Mono(nb_iter_tot).Resource_Matrix = Resource_Matrix;
    data_to_save_Mono(nb_iter_tot).Lag_time_Cons = Lag_time_Cons;
    data_to_save_Mono(nb_iter_tot).Mat_kappa_3 = Mat_kappa_3;
    data_to_save_Mono(nb_iter_tot).kappa_mat = kappa_mat;
    data_to_save_Mono(nb_iter_tot).R = R;
    data_to_save_Mono(nb_iter_tot).Death_Mat_init = Death_Mat_init;
    data_to_save_Mono(nb_iter_tot).Resource_Matrix_init = Resource_Matrix_init;
    data_to_save_Mono(nb_iter_tot).Lag_time_Cons_init = Lag_time_Cons_init;
    data_to_save_Mono(nb_iter_tot).Mat_kappa_3_init = Mat_kappa_3_init;
    data_to_save_Mono(nb_iter_tot).kappa_mat_init = kappa_mat_init;
    data_to_save_Mono(nb_iter_tot).Struct_Res_Mat = Struct_Res_Mat;
    data_to_save_Mono(nb_iter_tot).Struct_Death_Mat = Struct_Death_Mat;
    data_to_save_Mono(nb_iter_tot).death_rate = death_rate;
end

save(strcat(cd, '/Data/', 'data_to_save_SynCom21_Reduced'), 'data_to_save_Mono')