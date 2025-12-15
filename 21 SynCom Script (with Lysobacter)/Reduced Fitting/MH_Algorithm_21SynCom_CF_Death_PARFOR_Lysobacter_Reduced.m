%Script made by Adline Vouillamoz and Isaline Guex. Metropolis-Hasting
%algorithm to fit the co-culture interactions. The fitted parameters are the
%CF rates, the resources consumptions rates (initial value based on
%observed mono-cultures) and the self inhibition rates. Possibility to use
%a standard loop or a parallele loop (parfor)

%%% In this script we seek to fit only the interaction for Lysobacter,
%%% keeping the other parameters equal to SynCom20 fitting
clear;
close all;

%Save or Not
save_data = 1; %1 if save, 0 otherwise %Do not save for now
% Create a structure or cell array to hold the data
data_Reduced_dim = struct(); %To avoid overwriting, change the name of the save variable at the last line of the script. Name it accoridng to this structure name.
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')

%Parameters drawn from SnyCom20 fitting
Name_file = 'data_to_save_Inter_v2';%'data_to_save_with_corrections';%'data_to_save_5.mat';%'data_to_save.mat';%'Fixed_Res_12_Groups_v1'; %'Death_CF_Model_New_Mu_max_v10';%'Death_CF_Model_V17';%'Rand_Parameters';%
load(strcat('Data/', Name_file));

% parfor nb_iter_tot = 1:20   
for nb_iter_tot = 1:1

    %Loadind data
    Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
    Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');
    Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','142:163', 'Format','auto');
    Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 10, 'Range','1:22', 'Format','auto');%readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto');
    % Data_Evol_test = Data_Evol;
    Time_step = [0 12 22 24 38 70 72 168 240 504]; %[0 1 3 7 10 21]*24; %Measured time step in hours %CONTROL IT
    tspan = [0, max(Time_step)]; %Time interval in hours
    S = height(Data_Evol);

    %Turn 7 into 5
    yield_vect = table2array(Parameters_set(1:S,7)); %Change it according to the desired mu_max
    name = string(table2array(Parameters_set(1:S,1)));%string(table2array(table(2:22,3)));
    
    mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8)); %Log-normal Parameters_Senka_mu_max(:,7:8)
    mu_max_vect = mu_max_dist(:,1);
    Parameters_Senka_Lag_time = table2array(Parameters_Senka_Lag_time(:,7:8)); %table2array(Parameters_Senka_Lag_time(:,2:3));
    
    Data_Evol_temp = table2array(Data_Evol(:, 2:end));
    %Remove the second replicate for the fitting
    Data_Evol_test = Data_Evol_temp(:, ismember(mod(1:length(Data_Evol_temp(1,:)), 4), 2));
    Data_Evol_temp = Data_Evol_temp(:,~ismember(mod(1:length(Data_Evol_temp(1,:)), 4), 2));%Data_Evol_temp(:,~ismember(mod(1:length(Data_Evol_temp(1,:)), 8),[1, 2, 3]));%Data_Evol_temp(:,mod(1:length(Data_Evol_temp(1,:)),4) ~= 0);%80% of the data to train. Hold-out method.
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
    G_target = nb_time_step*4;   %about 4× more pseudo-points than observations
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
 
    %Load parameters
    data_to_save = data_to_save_Inter_v2;%data_to_save_with_corrections;%%data_to_save_5;%data_to_save_4;
    % Fitted parameters
    ind_sim = 3; %3, 5, (7), 9, (10), 11, 13, (15), (17), 18
    kappa_mat = data_to_save(ind_sim).kappa_mat;
    CrossFeed_Mat_Temp = data_to_save(ind_sim).CrossFeed_Mat_Temp;
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix;
    Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; 
    Threshold_CF = data_to_save(ind_sim).Threshold_CF;
    Threshold_death = data_to_save(ind_sim).Threshold_death; 
    Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
    R = data_to_save(ind_sim).R;
    nb_Res = length(R); %Number of resources (group of resources)
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp./kappa_mat(:,2);
    ind_Lyso = 6;

    death_rate = -[0;-0.000917642333192334;0;0;-0.000631406110652123; -0.000354932;-0.00122911113259887;...
        -0.00142098782776364;-0.00168390810416051;-0.000942517121033203;0;-0.000720301703810259;0;...
        0;-0.000192751969464189;-0.00254373051637326;-0.00163181578562513;0;-0.00132840048140091;-0.00185805721783005;-0.00623013508583033];
    death_rate_2 = [1e-05; 6.395e-04; 3.5730e-04; 7.1931e-04;  8.2640e-05; 3.54932e-04; 6.261921e-04;  4.6248e-05; 2.6131e-04; 3.4113e-04;...
    9.4029e-04; 9.9890e-05; 3.8416e-04; 1e-05; 6.0388e-06; 5.8124e-04; 2.1584e-04; 6.5872e-04; 2.6565e-04; 1.01239e-03; 1.30413e-03];
    death_rate = mean([death_rate, death_rate_2], 2);
    death_rate = max(death_rate, 7e-05);
    mu_max_vect_temp = exp(mu_max_vect + mu_max_dist(:,2).^2./2);
    yield = max(normrnd(yield_vect, 0.0), 0);
    kappa_mat = Increased_Mat(kappa_mat, ind_Lyso, S, 4, 0);
    kappa_mat(ind_Lyso,:) = [2.5e+05 mu_max_vect_temp(6) (mu_max_vect_temp(ind_Lyso)/yield(ind_Lyso) - mu_max_vect_temp(ind_Lyso)) 2.5e+05];
    Death_Mat_Temp = Increased_Mat(Death_Mat_Temp, ind_Lyso, S, S, 1); 
    Death_Mat_Temp(ind_Lyso, ind_Lyso) = 0.0214;
    CrossFeed_Mat_Temp = Increased_Mat(CrossFeed_Mat_Temp, ind_Lyso, S, S, 1);
    CrossFeed_Mat_Temp(ind_Lyso, :) = mean(CrossFeed_Mat_Temp([17 19], :));
    CrossFeed_Mat_Temp(:, ind_Lyso) = mean(CrossFeed_Mat_Temp(:,[17 19]), 2);
    CrossFeed_Mat_Temp(ind_Lyso, ind_Lyso) = 0;
    Lag_time_Cons = Increased_Mat(Lag_time_Cons, ind_Lyso, S, nb_Res, 0);
    Lag_time_Cons(ind_Lyso,:) = lognrnd(Parameters_Senka_Lag_time(ind_Lyso,1), Parameters_Senka_Lag_time(ind_Lyso,2));
    x = rand(1, nb_Res);
    rand_indice = x > 0.3;
    Lag_time_Rand = normrnd(50, 10, 1, nb_Res);
    Lag_time_Cons(ind_Lyso, rand_indice) = max(Lag_time_Cons(ind_Lyso, rand_indice), Lag_time_Rand(rand_indice));
    Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
    x = rand(1, nb_Res);
    rand_inidice = x > 0.24;
    Resource_Matrix = Increased_Mat(Resource_Matrix, ind_Lyso, S, nb_Res, 0); 
    Resource_Matrix = Resource_Matrix.*Var_Resource_Matrix;
    Resource_Matrix(ind_Lyso, rand_inidice) = lognrnd(zeros(1, sum(rand_inidice)), mu_max_dist(ind_Lyso,2).*ones(1, sum(rand_inidice)));
    Mat_kappa_3_Lyso = kappa_mat(:,3).*CrossFeed_Mat_Temp./kappa_mat(:,2);
    Threshold_CF = Increased_Mat(Threshold_CF, ind_Lyso, S, 1, 0);
    Threshold_CF(ind_Lyso) = 1.0e-03;
    Threshold_death = Increased_Mat(Threshold_death, ind_Lyso, S, 1, 0);
    Threshold_death(ind_Lyso) = 1.2e-04;
    
    liste_nrj_tot_temp = 1e15;
    
    %Setting for Matlab ODE solver
    opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
    
    %Initialization of the colors
    
    colors = distinguishable_colors(S);
    
    num_fig = 1;
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0;
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];
    
    %Fit one unique threshold value
    Threshold_Pred = max(normrnd(1e-04, 0, S, 1),0);%0?
    mat_y_0 = reshape(mat_y_0', 1, []);
    mat_y_0 = [mat_y_0 R];
    kappa_mat_init = kappa_mat; 
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp./mu_max_vect_temp;
    Mat_kappa_3_init = Mat_kappa_3;
    %Initialization predation matrix Lysobacter
    Pred_Mat_Lyso = zeros(S,S);
    Prey_num = [7 15 16 17];%[11 15 17 21];%[11 15 17];%[11 17 21];%[3 7 11 13 15 17 20 21];
    Pred_Mat_Lyso(ind_Lyso, Prey_num) = 1e-02; %Predation rate randomly chosen
    Binary_Pred = zeros(S); Binary_Pred(ind_Lyso, Prey_num) = 1;
    yield_Pred = 0.2;% 20% of yield for predation
    Lag_time_Pred = zeros(S, S);
    Lag_time_Pred(ind_Lyso, Prey_num) = 100;%72; %Assume the predation starts only when the primary resources are depleted
    Lag_time_Pred_init = Lag_time_Pred;
        
    beta = 0.5*1e2; 
    num_steps = 100;%500;
    nb_loops = 4; %Number of different loops (CF, death, resources, predation), 3 if predation loop is commented
         
    n = 1; r = 1; l = 1; p = 1; u = 1; h = 1; %Total energy
    max_tot_iter = 10;%15;%13;%18;%
    liste_nrj_CF = zeros(1, num_steps*max_tot_iter);
    Struct_CrossFeed_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
    Struct_CrossFeed_Mat{1} = CrossFeed_Mat_Temp; CrossFeed_Mat_init = CrossFeed_Mat_Temp;
    CrossFeed_Mat_Temp_cand = CrossFeed_Mat_Temp; %Initilization candidat
    Threshold_CF_cand = Threshold_CF; Threshold_CF_init = Threshold_CF;
    acceptance_ratio_CF = zeros(1, num_steps*max_tot_iter);
    liste_nrj_pred = zeros(1, num_steps*max_tot_iter);
    Struct_Death_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
    Struct_Death_Mat{1} = Death_Mat_Temp; Death_Mat_init = Death_Mat_Temp;
    Death_Mat_Temp_cand = Death_Mat_Temp; %Initilization candidat
    Threshold_death_cand = Threshold_death; Threshold_death_init = Threshold_death_cand;
    acceptance_ratio_pred = zeros(1, num_steps*max_tot_iter);
    liste_nrj_res = zeros(1, num_steps*max_tot_iter);
    Struct_Res_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
    Struct_Res_Mat{1} = Resource_Matrix; Resource_Matrix_init = Resource_Matrix;
    Resource_Matrix_cand = Resource_Matrix; %Initilization candidat
    acceptance_ratio_res = zeros(1, num_steps*max_tot_iter);
    %%%
    Struct_Pred_Mat_Lyso = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
    Struct_Pred_Mat_Lyso{1} = Pred_Mat_Lyso; Pred_Mat_Lyso_init = Pred_Mat_Lyso;
    Pred_Mat_Lyso_cand = Pred_Mat_Lyso; %Initilization candidat
    %%%
    liste_nrj_tot = zeros(1, nb_loops*num_steps*max_tot_iter);
    acceptance_ratio_tot = zeros(1, nb_loops*num_steps*max_tot_iter);
    nb_Res_new = nb_Res;
    Lag_time_Cons_cand = Lag_time_Cons; Lag_time_Cons_init = Lag_time_Cons;
    R_new = R; R_init = R; mat_y_0_new = mat_y_0; 
    liste_nrj_nb_res = zeros(1, num_steps*max_tot_iter);
    acceptance_ratio_nb_res = zeros(1, num_steps*max_tot_iter);
    Evol_nb_Res = zeros(1, num_steps*max_tot_iter);
    
    var_theta = 0.1;
    acceptance_ratio_temp = 10000;
    ratio_increase_loop = 1;
    N = nb_loops*num_steps*max_tot_iter*round(ratio_increase_loop^max_tot_iter);
    max_val = [10./mu_max_vect_temp 10*ones(S,1) 10./mu_max_vect_temp];
    Hill_CF_coeff = 10; Hill_Pred_coeff = 10;
    Hill_CF_coeff_cand = Hill_CF_coeff; Hill_Pred_coeff_cand = Hill_Pred_coeff;
    delta_dist = 1;
    LT_CF_Death = zeros(S, S);
    Cov_LT_eye = eye(S*S);
    Weight_species =  ones(S,8);
    weight_day_seven = 1;
    CF_Trace_Plot = []; Death_Trace_Plot = []; Res_Trace_Plot = [];
    Cov_Mat_CF = 0.01*eye(S*S); Cov_Mat_Death = eye(S); Cov_Mat_Res = 0.01*eye(S*nb_Res);
    lag = 5;
    [param_chain_1, param_chain_2, param_chain_4] = deal(zeros(S, S, num_steps*max_tot_iter));
    param_chain_3 = zeros(S, nb_Res, num_steps*max_tot_iter);
    while p < N
    
        %%%%%%%%%% Cross-feeding loop %%%%%%%% 
    
        theta_T = 0.01*var_theta*mean(Threshold_CF)*eye(S); %Remove eye(S)
        Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp_cand./mu_max_vect_temp;%CrossFeed_Mat_Temp_cand./yield' - CrossFeed_Mat_Temp_cand;
        Cov_Mat = 0.1*eye(S*S);%Cov_Mat_CF;%
        W_n_Cov_Mat = normrnd(zeros(S*S,1), ones(S*S,1)); W_n_LT = normrnd(zeros(S*S,1), ones(S*S,1)); 
        accept_temp = zeros(1, num_steps);
        Temp = 500;
        for k = 1:num_steps

            for j = 1:lag
                sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat_Temp_cand, Mat_kappa_3, Resource_Matrix,...
                    Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate,...
                    Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); %Multiple resource groups
                z_temp = deval(sol, Time_step);
                X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
                P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
                X = X + P;
                nu = k^(-3/4);
                StackPlotTot = X./sum(X);

                % === GRADIENT MATCHING PENALTY ===
                lambda_grad = 1e-1;  % tune between 1e-2 and 1
                
                S = size(Measured_Abund, 1);
                f_theta = zeros(S, numel(t_grid));
                
                for g = 1:numel(t_grid)
                    % Evaluate state from ODE solution
                    y_curr = deval(sol, t_grid(g));  % same as in your code below
                    % Compute derivative using your existing ODE RHS
                    dz_curr = fun_CF_Death_Lyso(t_grid(g), y_curr, kappa_mat, CrossFeed_Mat_Temp_cand, Mat_kappa_3, Resource_Matrix, ...
                        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate, Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, ...
                        Lag_time_Pred, nb_Res, 10, 10);
                    % Extract only the dx part (species biomass derivative)
                    f_theta(:, g) = dz_curr(1:S);
                end
                
                % Compute gradient mismatch
                diff_grad = (dx_hat - f_theta).^2 ./ (var_dx + 1e-8);
                E_grad = lambda_grad * sum(diff_grad, 'all');
                %%%
    
                CrossFeed_Mat_Temp = CrossFeed_Mat_Temp./(mu_max_vect_temp);
                CrossFeed_Mat_Temp_cand = CrossFeed_Mat_Temp_cand./(mu_max_vect_temp);
        
                [accept_temp, W_n_LT, W_n_Cov_Mat, CrossFeed_Mat_Temp_cand, CrossFeed_Mat_Temp, Struct_CrossFeed_Mat, liste_nrj_tot, liste_nrj_CF, acceptance_ratio_tot, acceptance_ratio_CF, p, n, nu, Cov_Mat, param_chain_1] = ...
                    fun_MH_Candidate_Rob(accept_temp, k, W_n_LT, W_n_Cov_Mat, CrossFeed_Mat_Temp_cand, CrossFeed_Mat_Temp, Struct_CrossFeed_Mat, p, n, liste_nrj_tot, liste_nrj_CF,...
                    acceptance_ratio_tot, acceptance_ratio_CF, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, nu, Cov_Mat, Cov_LT_eye, max_val(:,1), nb_rep, LT_CF_Death, LT_CF_Death, Weight_species, weight_day_seven, Temp, param_chain_1, E_grad);

                %%%
                temp = CrossFeed_Mat_Temp;
                CrossFeed_Mat_Temp_cand(1:1+size(CrossFeed_Mat_Temp_cand,1):end) = 0;
                temp(ind_Lyso,:) = CrossFeed_Mat_Temp_cand(ind_Lyso,:);
                temp(:,ind_Lyso) = CrossFeed_Mat_Temp_cand(:,ind_Lyso);
                CrossFeed_Mat_Temp_cand = temp;
                %%%
        
                CF_Trace_Plot(:, :, n) = CrossFeed_Mat_Temp;  
        
                CrossFeed_Mat_Temp_cand(CrossFeed_Mat_Temp_cand < 1e-06) = 0;

                CrossFeed_Mat_Temp = CrossFeed_Mat_Temp.*(mu_max_vect_temp);
                CrossFeed_Mat_Temp_cand = CrossFeed_Mat_Temp_cand.*(mu_max_vect_temp);
                Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp_cand./mu_max_vect_temp;
            end
            if mod(p, 50) == 0 && std(liste_nrj_tot(p - 49:p)) < 1e-8
                Temp = Temp*0.5;
            end

            Struct_CrossFeed_Mat{n + 1} = CrossFeed_Mat_Temp;
            n = n + 1;
            p = p + 1;
        
    
            if mod(p, num_steps) == 0
                disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
            end          
        end
        disp(sum(accept_temp)/num_steps)
        Cov_Mat_CF = Cov_Mat;

        burn_in = round(0.3*size(param_chain_1, 2));
        param_chain_post = param_chain_1(:, :, burn_in:n);
        
        % Compute mean, std, and CV for each parameter
        param_mean = mean(param_chain_post, 2);
        param_std  = std(param_chain_post, 0, 2);
        CV_1 = param_std./abs(param_mean);
        %%%%%%%%%% Death loop %%%%%%%%
        
        theta = var_theta*mean(mean(Death_Mat_Temp_cand)); 
        theta_T_death = 0.01*var_theta*mean(Threshold_death)*eye(S);
        Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp./mu_max_vect_temp;
        Cov_Mat = 0.1*eye(S);%Cov_Mat_Death;%
        W_n_thresh_death = normrnd(zeros(S,1), ones(S,1)); W_n_Cov_Mat_death = normrnd(zeros(S,1), ones(S,1)); 
        accept_temp = zeros(1, num_steps);
        Temp = 500;
        for k = 1:num_steps
    
            for j = 1:lag
                sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix,...
                    Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp_cand, death_rate,...
                    Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
                z_temp = deval(sol, Time_step);
                X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
                P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
                X = X + P;
                nu = k^(-3/4); 
                StackPlotTot = X./sum(X);

                % === GRADIENT MATCHING PENALTY ===
                lambda_grad = 1e-1;  % tune between 1e-2 and 1
                
                S = size(Measured_Abund, 1);
                f_theta = zeros(S, numel(t_grid));
                
                for g = 1:numel(t_grid)
                    % Evaluate state from ODE solution
                    y_curr = deval(sol, t_grid(g));  % same as in your code below
                    % Compute derivative using your existing ODE RHS
                    dz_curr = fun_CF_Death_Lyso(t_grid(g), y_curr, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix, ...
                        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp_cand, death_rate, Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, ...
                        Lag_time_Pred, nb_Res, 10, 10);
                    % Extract only the dx part (species biomass derivative)
                    f_theta(:, g) = dz_curr(1:S);
                end
                
                % Compute gradient mismatch
                diff_grad = (dx_hat - f_theta).^2 ./ (var_dx + 1e-8);
                E_grad = lambda_grad * sum(diff_grad, 'all');
                %%%
        
                [accept_temp, W_n_Cov_Mat_death, Death_Mat_Temp_cand, Death_Mat_Temp, Struct_Death_Mat, liste_nrj_tot, liste_nrj_pred, acceptance_ratio_tot, acceptance_ratio_pred, p, r, nu, Cov_Mat, param_chain_2] = ...
                    fun_MH_Candidate_Rob_Death(accept_temp, k, W_n_Cov_Mat_death, Death_Mat_Temp_cand, Death_Mat_Temp, Struct_Death_Mat, p, r, liste_nrj_tot, liste_nrj_pred,...
                    acceptance_ratio_tot, acceptance_ratio_pred, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, nu, Cov_Mat, max_val(:,2), nb_rep, Weight_species, weight_day_seven, Temp, param_chain_2, E_grad);
    
                temp = Death_Mat_Temp;
                temp(ind_Lyso,:) = Death_Mat_Temp_cand(ind_Lyso,:);
                temp(ind_Lyso, ind_Lyso) = 0;
                Death_Mat_Temp_cand = temp;

                Death_Trace_Plot(:,:, r) = Death_Mat_Temp;
                ind = eye(S,S);
                Death_Mat_Temp_cand(~ind) = 0; %Changement pour diagonale zeros
            end
            if mod(p, 50) == 0 && std(liste_nrj_tot(p - 49:p)) < 1e-8
                Temp = Temp*0.5;
            end

            Struct_Death_Mat{r + 1} = Death_Mat_Temp;
            r = r + 1;
            p = p + 1;
    
            if mod(p, num_steps) == 0
                disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
            end
        end
        disp(sum(accept_temp)/num_steps)
        Cov_Mat_Death = Cov_Mat;
    
        %%%%%%%%%% Resource matrix loop %%%%%%%%
    
        Cov_Mat = 0.1*eye(S*nb_Res);
        Cov_Mat_LT = 0.01*eye(S*nb_Res);
        W_n_LT_Res = normrnd(zeros(S*nb_Res,1), ones(S*nb_Res,1)); W_n_Cov_Mat_Res = normrnd(zeros(S*nb_Res,1), ones(S*nb_Res,1));
        accept_temp = zeros(1, num_steps);
        Temp = 500;
        for k = 1:(num_steps)

            for j = 1:(lag)
                sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix_cand,...
                    Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate,...
                    Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons_cand, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
                z_temp = deval(sol, Time_step);
                X = z_temp(mod(1:S*3,3) == 1, :); %Species'biomass
                P = z_temp(mod(1:S*3,3) == 2, :); %Complexes'biomass
                X = X + P;
                nu = k^(-3/4); 
                StackPlotTot = X./sum(X);

                % === GRADIENT MATCHING PENALTY ===
                lambda_grad = 1e-1;  % tune between 1e-2 and 1
                
                S = size(Measured_Abund, 1);
                f_theta = zeros(S, numel(t_grid));
                
                for g = 1:numel(t_grid)
                    % Evaluate state from ODE solution
                    y_curr = deval(sol, t_grid(g));  % same as in your code below
                    % Compute derivative using your existing ODE RHS
                    dz_curr = fun_CF_Death_Lyso(t_grid(g), y_curr, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix_cand, ...
                        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate, Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons_cand, ...
                        Lag_time_Pred, nb_Res, 10, 10);
                    % Extract only the dx part (species biomass derivative)
                    f_theta(:, g) = dz_curr(1:S);
                end
                
                % Compute gradient mismatch
                diff_grad = (dx_hat - f_theta).^2 ./ (var_dx + 1e-8);
                E_grad = lambda_grad * sum(diff_grad, 'all');
                %%%

                [accept_temp, W_n_LT_Res, W_n_Cov_Mat_Res, Resource_Matrix_cand, Resource_Matrix, Struct_Res_Mat, liste_nrj_tot, liste_nrj_res, acceptance_ratio_tot, acceptance_ratio_res, p, l, nu, Cov_Mat, param_chain_3, Cov_Mat_LT, Lag_time_Cons_cand, Lag_time_Cons] = ...
                    fun_MH_Candidate_Rob(accept_temp, k, W_n_LT_Res, W_n_Cov_Mat_Res, Resource_Matrix_cand, Resource_Matrix, Struct_Res_Mat, p, l, liste_nrj_tot, liste_nrj_res,...
                    acceptance_ratio_tot, acceptance_ratio_res, S, nb_Res, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, nu, Cov_Mat, Cov_Mat_LT, max_val(:,3), nb_rep, Lag_time_Cons_cand, Lag_time_Cons, Weight_species, weight_day_seven, Temp, param_chain_3, E_grad);
                        
                temp = Resource_Matrix;
                temp(ind_Lyso,:) = Resource_Matrix_cand(ind_Lyso,:);
                Resource_Matrix_cand = temp;
            end
            if mod(p, 50) == 0 && std(liste_nrj_tot(p - 49:p)) < 1e-8
                Temp = Temp*0.5;
            end

            Struct_Res_Mat{r + 1} = Resource_Matrix;
            l = l + 1;
            p = p + 1;


            if mod(p, num_steps) == 0
                disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
            end
        end
        disp(sum(accept_temp)/(num_steps))
        Cov_Mat_Res = Cov_Mat;

        %%%%%%%%%% Predation matrix loop %%%%%%%%
        %If kept modify the threshold too

        Cov_Mat = 0.1*eye(S*S);
        Cov_Mat_LT = 0.01*eye(S*S);
        W_n_LT_Pred = normrnd(zeros(S*S,1), ones(S*S,1)); W_n_Cov_Mat_Pred = normrnd(zeros(S*S,1), ones(S*S,1));
        accept_temp = zeros(1, num_steps);
        Temp = 500;
        for k = 1:(num_steps)

            for j = 1:(lag)
                sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix,...
                    Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate,...
                    Pred_Mat_Lyso_cand, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
                z_temp = deval(sol, Time_step);
                X = z_temp(mod(1:S*3,3) == 1, :); %Species'biomass
                P = z_temp(mod(1:S*3,3) == 2, :); %Complexes'biomass
                X = X + P;
                nu = k^(-3/4); 
                StackPlotTot = X./sum(X);

                % === GRADIENT MATCHING PENALTY ===
                lambda_grad = 1e-1;  %tune between 1e-2 and 1
                
                S = size(Measured_Abund, 1);
                f_theta = zeros(S, numel(t_grid));
                
                for g = 1:numel(t_grid)
                    % Evaluate state from ODE solution
                    y_curr = deval(sol, t_grid(g));  % same as in your code below
                    % Compute derivative using your existing ODE RHS
                    dz_curr = fun_CF_Death_Lyso(t_grid(g), y_curr, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix, ...
                        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate, Pred_Mat_Lyso_cand, yield_Pred, S, Lag_time_Cons, ...
                        Lag_time_Pred, nb_Res, 10, 10);
                    % Extract only the dx part (species biomass derivative)
                    f_theta(:, g) = dz_curr(1:S);
                end
                
                % Compute gradient mismatch
                diff_grad = (dx_hat - f_theta).^2 ./ (var_dx + 1e-8);
                E_grad = lambda_grad * sum(diff_grad, 'all');
                %%%

                [accept_temp, W_n_LT_Pred, W_n_Cov_Mat_Pred, Pred_Mat_Lyso_cand, Pred_Mat_Lyso, Struct_Res_Mat, liste_nrj_tot, liste_nrj_res, acceptance_ratio_tot, acceptance_ratio_res, p, l, nu, Cov_Mat, param_chain_4] = ...
                    fun_MH_Candidate_Rob(accept_temp, k, W_n_LT_Pred, W_n_Cov_Mat_Pred, Pred_Mat_Lyso_cand, Pred_Mat_Lyso, Struct_Res_Mat, p, l, liste_nrj_tot, liste_nrj_res,...
                    acceptance_ratio_tot, acceptance_ratio_res, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, nu, Cov_Mat, Cov_Mat_LT, max_val(:,3), nb_rep, LT_CF_Death, LT_CF_Death, Weight_species, weight_day_seven, Temp, param_chain_4, E_grad);

                Pred_Mat_Lyso_cand = Binary_Pred.*Pred_Mat_Lyso_cand; %multiplication by to get a 0-1 matrix. 1 at imposed prey 0 otherwise. We don't seek the preys, we assume there are known. To discuss but otherwise I'm afraid there would have too many possibilites (already true)
            end
            if mod(p, 50) == 0 && std(liste_nrj_tot(p - 49:p)) < 1e-8
                Temp = Temp*0.5;
            end

            Struct_Pred_Mat_Lyso{r + 1} = Pred_Mat_Lyso;
            h = h + 1;
            p = p + 1;


            if mod(p, num_steps) == 0
                disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
            end
        end
        disp(sum(accept_temp)/(num_steps))
        Cov_Mat_Res = Cov_Mat;
        %%%
    
        Resource_Matrix_cand = Resource_Matrix; nb_Res_new = nb_Res; Lag_time_Cons_cand = Lag_time_Cons; R_new = R; mat_y_0_new = mat_y_0; %Re-initialization if the last candidate is not accepted. Otherwise issues with the resource loop.
        acceptance_ratio_temp = sum(acceptance_ratio_tot)/p;
        liste_nrj_tot_temp = liste_nrj_tot(p - 1);
        delta_dist = 1;
        num_steps = round(ratio_increase_loop*num_steps); 
        weight_day_seven = ratio_increase_loop*weight_day_seven;
        ratio_obs_sim = 1 - abs(X - mean(Measured_Abund,3))./(abs(X) + abs(mean(Measured_Abund,3))); %Similaritiy intex used as weights for variances
        Weight_species = mean(ratio_obs_sim, 2);
    end 
    acceptance_ratio_temp
    liste_nrj_tot_temp
   
    
    %Data to test on half of the remaining data
    % Time_step = [0 1 3 7 10 21]*24; %Measured time step in hours. If
    % another time step on the test set
    tspan = [0, max(Time_step)]; %Time interval in hours
    Data_Evol_temp = Data_Evol_test; %table2array(Data_Evol_test(:, 2:end));
    nb_time_step = length(Time_step);
    nb_obs = length(Data_Evol_temp(1,:));
    nb_rep = nb_obs/nb_time_step;
    Measured_Abund = zeros(length(mu_max_vect), nb_time_step, nb_rep); %Number species, number times, number replicates.
    for i = 1:nb_rep
        Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
    end
    StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
    Threshold_Surviving = 1e-10;
    nb_Surv_Obs = sum(mean(Measured_Abund,3) > Threshold_Surviving);
    mean_y_0 = mean(Measured_Abund(:,1,:), 3);%table2array(Data_Evol(1:20, 2));
    mat_y_0 = [mean_y_0 zeros(S,1) zeros(S,1)];
    mat_y_0 = reshape(mat_y_0', 1, []);
    mat_y_0 = [mat_y_0 R];
    
    sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix,...
        Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate,...
        Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
    z_temp = deval(sol, Time_step);
    X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
    P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
    X = X + P;
    
            
    %%% Figure generation %%%
    
    acceptance_ratio_fin = sum(acceptance_ratio_tot)/(nb_loops*num_steps*max_tot_iter);
    ratio_fin_abund = sum(X(:,end))/sum(Measured_Abund(:,end));
    
    % figure(num_fig) %figure 1 : Simulated absolute abundances
    for j = 1:S
        plot(Time_step, X(j,:), '-o', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
        hold on
    end
    axis([0 500 0 3.5e-04])
    %legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
    
    hold on
    %num_fig = num_fig + 1;
    figure(num_fig) %figure 2 : Observed absolute abundances
    for j = 1:S
        plot(Time_step, mean(Measured_Abund(j,:,:), 3), '--*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
        hold on
    end
    axis([0 500 0 3.5e-04])
    legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
    
    num_fig = num_fig + 1; 
    z_temp = z_temp(1:(end-nb_Res), end);
    z_temp = reshape(z_temp',3, S);
    z_temp = z_temp';
    z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
    StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
    [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
    [Shann_obs, Simp_obs] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));
    rand_prop = max(normrnd(1, 0.3, S, 1), 0);
    mat_y_0 = (z_temp(:,1:2)/10).*rand_prop; %Take the 10% of the system
    mat_y_0 = [mat_y_0 zeros(S,1)];
    
    StackPlotTot = X./sum(X);
    
    z_fin_sim = z(1:S,end); %Absolute stationary abundances
    z_fin_obs = mean(Measured_Abund(1:S, end, :), 3); %Absolute stationary abundances
    nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
    sum(z_fin_sim)/sum(z_fin_obs)
    acceptance_ratio_temp
    liste_nrj_tot_temp
    
    figure(num_fig); 
    bar(StackPlotTot', 'stacked');
    axis([0 11.5 0 1])
    legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
    title('Stacked all replicates')
    num_fig = num_fig + 1;
    
    figure(num_fig); 
    bar(mean(StackPlot_Meas,3)', 'stacked');
    axis([0 11.5 0 1])
    legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
    title('Stacked observed')
    num_fig = num_fig + 1;
    
    [diff_vect, I, ~] = Sort_diff(StackPlotTot(:,end),mean(StackPlot_Meas(:,end,:),3));
    name_order = name(I);
    figure(num_fig);
    stem(diff_vect);
    xtickangle(90)
    set(gca,'xtick',1:21,'xticklabel',name_order)
    ylabel('Sim - Obs')
    
    num_fig = num_fig + 1;
    
    figure(num_fig);
    for j = 1:S
        scatter(mean(StackPlot_Meas(j, 2, :), 3), StackPlotTot(j,2), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors(j,:), 'MarkerFaceColor', colors(S+1-j,:));%col(S+1-j,:))      
        hold on
    end
    axis([0 0.5 0 0.5]);
    axis square
    reflin = refline(1,0);
    reflin.Color = 'r';
    xlabel('Experiment'); 
    ylabel('Simulation');
    legend(name);
    title('Scatter Tot');
    
    num_fig = num_fig + 1;
    
    figure(num_fig);
    for j = 1:S
        scatter(z_fin_obs(j), z_fin_sim(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors(j,:), 'MarkerFaceColor', colors(S+1-j,:));%col(S+1-j,:))      
        hold on
    end
    axis([0 3*10^(-4) 0 3*10^(-4)]);
    reflin = refline(1,0);
    axis square
    reflin.Color = 'r';
    xlabel('Experiment'); 
    ylabel('Simulation');
    legend(name);
    title('Scatter absolute abundance');
    num_fig = num_fig + 1;
    
    figure(num_fig);
    plot(1:length(Time_step), Simp_obs, 'b--o')
    hold on
    plot(1:length(Time_step), Simp_sim, 'r--o')
    num_fig = num_fig + 1;
    
    figure(num_fig);
    plot(1:length(Time_step), nb_Surv_Obs, 'b--o')
    hold on
    plot(1:length(Time_step), nb_Surv_Sim, 'r--o')
    num_fig = num_fig + 1;
    
    num_fig = num_fig+1;
    figure(num_fig);
    plot(1:length(liste_nrj_tot), liste_nrj_tot)

    data_Reduced_dim(nb_iter_tot).CrossFeed_Mat_Temp = CrossFeed_Mat_Temp;   
    data_Reduced_dim(nb_iter_tot).Pred_Mat_Lyso = Pred_Mat_Lyso;
    data_Reduced_dim(nb_iter_tot).Threshold_Pred = Threshold_Pred;
    data_Reduced_dim(nb_iter_tot).Threshold_CF = Threshold_CF;
    data_Reduced_dim(nb_iter_tot).Threshold_death = Threshold_death;
    data_Reduced_dim(nb_iter_tot).Death_Mat_Temp = Death_Mat_Temp;
    data_Reduced_dim(nb_iter_tot).Resource_Matrix = Resource_Matrix;
    data_Reduced_dim(nb_iter_tot).Lag_time_Pred = Lag_time_Pred;
    data_Reduced_dim(nb_iter_tot).Lag_time_Cons = Lag_time_Cons;
    data_Reduced_dim(nb_iter_tot).Mat_kappa_3 = Mat_kappa_3;
    data_Reduced_dim(nb_iter_tot).kappa_mat = kappa_mat;
    data_Reduced_dim(nb_iter_tot).R = R;
    data_Reduced_dim(nb_iter_tot).CrossFeed_Mat_init = CrossFeed_Mat_init;
    data_Reduced_dim(nb_iter_tot).Threshold_CF_init = Threshold_CF_init;
    data_Reduced_dim(nb_iter_tot).Pred_Mat_Lyso_init = Pred_Mat_Lyso_init;
    data_Reduced_dim(nb_iter_tot).Threshold_death_init = Threshold_death_init;
    data_Reduced_dim(nb_iter_tot).Death_Mat_init = Death_Mat_init;
    data_Reduced_dim(nb_iter_tot).Resource_Matrix_init = Resource_Matrix_init;
    data_Reduced_dim(nb_iter_tot).Lag_time_Pred_init = Lag_time_Pred_init;
    data_Reduced_dim(nb_iter_tot).Lag_time_Cons_init = Lag_time_Cons_init;
    data_Reduced_dim(nb_iter_tot).Mat_kappa_3_init = Mat_kappa_3_init;
    data_Reduced_dim(nb_iter_tot).kappa_mat_init = kappa_mat_init;
    data_Reduced_dim(nb_iter_tot).R_init = R_init;
    data_Reduced_dim(nb_iter_tot).Struct_Res_Mat = Struct_Res_Mat;
    data_Reduced_dim(nb_iter_tot).Struct_Death_Mat = Struct_Death_Mat;
    data_Reduced_dim(nb_iter_tot).Struct_CrossFeed_Mat = Struct_CrossFeed_Mat;
    data_Reduced_dim(nb_iter_tot).Struct_Pred_Mat_Lyso = Struct_Pred_Mat_Lyso;
    data_Reduced_dim(nb_iter_tot).Cov_Mat_CF = Cov_Mat_CF;
    data_Reduced_dim(nb_iter_tot).Cov_Mat_Death = Cov_Mat_Death;
    data_Reduced_dim(nb_iter_tot).Cov_Mat_Res = Cov_Mat_Res;
    data_Reduced_dim(nb_iter_tot).Cov_Mat_LT = Cov_Mat_LT;
    data_Reduced_dim(nb_iter_tot).death_rate = death_rate;
end

save(strcat(cd, '/Data/', 'data_Reduced_dim'), 'data_Reduced_dim')