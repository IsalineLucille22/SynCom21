%Script made by Adline Vouillamoz and Isaline Guex. Metropolis-Hasting
%algorithm to fit the co-culture interactions. The fitted parameters are the
%CF rates, the resources consumptions rates (initial value based on
%observed mono-cultures) and the self inhibition rates. Possibility to use
%a standard loop or a parallele loop (parfor)

%%% To do; fit parameters for the SynCom21 with Senka's data. Then test it
%%% (in Main_Script_SynCom21.m) on Clara's data. MAKE A LOOP FOR THE
%%% PREDATION MATRIX?
clear;
close all;

%Save or Not
save_data = 1; %1 if save, 0 otherwise %Do not save for now
% Create a structure or cell array to hold the data
data_to_save_SynCom21_Soil_with_pred_OneRes_Parfor = struct(); %To avoid overwriting, change the name of the save variable at the last line of the script. Name it accoridng to this structure name.
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')
% parfor nb_iter_tot = 1:4   
for nb_iter_tot = 1:1

    %Loadind data
    Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
    Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');
    Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','142:163', 'Format','auto');
    Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto');
    Data_Evol_test = Data_Evol;
    Time_step = [0 1 3 7 10 21]*24; %Measured time step in hours %CONTROL IT
    tspan = [0, max(Time_step)]; %Time interval in hours
    S = height(Data_Evol);

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
    Data_Evol_temp = Data_Evol_temp(:,~ismember(mod(1:length(Data_Evol_temp(1,:)), 4), 2));%Data_Evol_temp(:,~ismember(mod(1:length(Data_Evol_temp(1,:)), 8),[1, 2, 3]));%Data_Evol_temp(:,mod(1:length(Data_Evol_temp(1,:)),4) ~= 0);%80% of the data to train. Hold-out method.
    % Data_Evol_temp(:, 1:4) = correction_table.*Data_Evol_temp(:, 1:4);
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
    
    liste_nrj_tot_temp = 1e15;
    %Initialization of the model parameters fixed for all replicates 
    t_0 = 0; %Time 0
    nb_Res = 1;%12; %Number of resources (group of resources)
    nb_Res_max = 1;%12;
    CrossFeed_Mat = lognrnd((log(0.8) + mu_max_dist(:,1)).*ones(S,S), 0*mu_max_dist(:,2).*ones(S,S));%.val val% of zeros. ConsumerxProducer
    Lag_time_Pred = zeros(S, S);
    death_rate = -[0;-0.000917642333192334;0;0;-0.000631406110652123; -0.000354932;-0.00122911113259887;...
        -0.00142098782776364;-0.00168390810416051;-0.000942517121033203;0;-0.000720301703810259;0;...
        0;-0.000192751969464189;-0.00254373051637326;-0.00163181578562513;0;-0.00132840048140091;-0.00185805721783005;-0.00623013508583033];
    death_rate_2 = [1e-05; 6.395e-04; 3.5730e-04; 7.1931e-04;  8.2640e-05; 3.54932e-04; 6.261921e-04;  4.6248e-05; 2.6131e-04; 3.4113e-04;...
    9.4029e-04; 9.9890e-05; 3.8416e-04; 1e-05; 6.0388e-06; 5.8124e-04; 2.1584e-04; 6.5872e-04; 2.6565e-04; 1.01239e-03; 1.30413e-03];
    death_rate = mean([death_rate, death_rate_2], 2);
    CrossFeed_Mat_Temp = zeros(S,S); %Temporal predation matrix
    rand_indice = rand(S,S) > 0.8; %Percentage of zeros. Larger is the value smaller is the number of interaction
    CrossFeed_Mat_Temp(rand_indice) = CrossFeed_Mat(rand_indice); %Put some element to 0
    Death_Mat_Temp = zeros(S,S);
    Death_Mat_Temp(1:1+size(Death_Mat_Temp,1):end) = 0.6*death_rate;
    death_rate = max(death_rate, 7e-05);
    CrossFeed_Mat_Temp(1:1+size(CrossFeed_Mat_Temp,1):end) = 0;
    
    % Create a resource matrix
    %Load parameters
    R_0 = 0.5;
    x = rand(S,nb_Res);
    Res_Percentage_old = Res_Percentage;
    Res_Percentage = min(max(Res_Percentage, 0.24), 0.24);
    rand_indice = x > Res_Percentage;
    Resource_Matrix = ones(S,nb_Res);
    Resource_Matrix(rand_indice) = 0;
    Resource_Matrix = Resource_Matrix.*lognrnd(-mu_max_dist(:,2).^2./2.*ones(S, nb_Res), mu_max_dist(:,2).*ones(S, nb_Res));%Resource_Matrix.*max(normrnd(ones(S, nb_Res), 0.05.*ones(S, nb_Res)), 0.01);%Resource_Matrix.*max(normrnd(ones(S, nb_Res), mu_max_dist(:,2)./mu_max_dist(:,1).*ones(S, nb_Res)), 0.01);%Centered on mu_max because multiplied mby mu_max so depend on the measured mu_max.
    R = R_0*15*10^(-3);%.*[0.6990 0.0581 0.0029 0.1477 4.0011e-05 0.0061 0.0233 0.0288 0.0008 0.0168 0.0070 0.0094]; %based on the observed values
    Lag_time_Cons = lognrnd(Parameters_Senka_Lag_time(:,1).*ones(S,nb_Res), Parameters_Senka_Lag_time(:,2).*ones(S,nb_Res)); %max(normrnd(Parameters_Senka_Lag_time(:,1).*ones(S,nb_Res), Parameters_Senka_Lag_time(:,2).*ones(S,nb_Res)), 0); %Lag time for resource consumption
    x = rand(S,nb_Res);
    rand_indice = x > 0.3; %Res_Percentage_Nb_Cons;
    Lag_time_Rand = normrnd(50, 10, S, nb_Res);%lognrnd(log(50), 1, S, nb_Res);
    Lag_time_Cons(rand_indice) = max(Lag_time_Cons(rand_indice), Lag_time_Rand(rand_indice));
    Res_Percentage = Res_Percentage_old;
    
    %Setting for Matlab ODE solver
    opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
    
    %Initialization of the colors
    
    colors = distinguishable_colors(S);
    
    num_fig = 1;
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0;
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];
    
    %Fit one unique threshold value
    Threshold_CF = max(normrnd(1.0e-03, 0*0.1e-05, S, 1),0); 
    Threshold_death = max(normrnd(1.2e-04, 0*0.1e-05, S, 1),0);
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
    %Initialization predation matrix Lysobacter
    Pred_Mat_Lyso = zeros(S,S);
    Prey_num = [7 15 16 17];%[11 15 17 21];%[3 7 11 13 15 17 20 21];
    ind_Lyso = 6;
    Pred_Mat_Lyso(ind_Lyso, Prey_num) = 1e-02; %Predation rate randomly chosen
    Pred_Mat_Lyso(ind_Lyso, 21) = 1;
    Binary_Pred = zeros(S); Binary_Pred(ind_Lyso, Prey_num) = 1;
    yield_Pred = 0.2;% 20% of yield for predation
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
    theta_res = var_theta*mean(mu_max_dist(:,2));
    delta_dist_vect = [];
    LT_CF_Death = zeros(S, S);
    Cov_Lt_eye = eye(S*S);
    Weight_species =  ones(S,8);
    weight_day_seven = 1;
    CF_Trace_Plot = []; Death_Trace_Plot = []; Res_Trace_Plot = [];
    Cov_Mat_CF = 0.01*eye(S*S); Cov_Mat_Death = eye(S); Cov_Mat_Res = 0.01*eye(S*nb_Res);
    lag = 5;
    while p < N
    
        %%%%%%%%%% Cross-feeding loop %%%%%%%% 
    
        theta_T = 0.01*var_theta*mean(Threshold_CF)*eye(S); %Remove eye(S)
        Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp_cand./mu_max_vect_temp;%CrossFeed_Mat_Temp_cand./yield' - CrossFeed_Mat_Temp_cand;
        Cov_Mat = 0.1*eye(S*S);%Cov_Mat_CF;%
        W_n_Cov_Mat = normrnd(zeros(S*S,1), ones(S*S,1)); W_n_LT = normrnd(zeros(S*S,1), ones(S*S,1)); W_n_thresh = normrnd(zeros(S,1), ones(S,1));  
        accept_temp = zeros(1, num_steps);
        Temp = 500;%1;
        for k = 1:num_steps

            for j = 1:lag
                sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat_Temp_cand, Mat_kappa_3, Resource_Matrix,...
                    Threshold_CF_cand, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate,...
                    Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); %Multiple resource groups
                z_temp = deval(sol, Time_step);
                X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
                P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
                X = X + P;
                nu = k^(-3/4);
                StackPlotTot = X./sum(X);
    
                CrossFeed_Mat_Temp = CrossFeed_Mat_Temp./(mu_max_vect_temp);
                CrossFeed_Mat_Temp_cand = CrossFeed_Mat_Temp_cand./(mu_max_vect_temp);
        
                [accept_temp, W_n_LT, W_n_thresh, W_n_Cov_Mat, CrossFeed_Mat_Temp_cand, CrossFeed_Mat_Temp, Threshold_CF_cand, Threshold_CF, Struct_CrossFeed_Mat, liste_nrj_tot, liste_nrj_CF, acceptance_ratio_tot, acceptance_ratio_CF, p, n, theta_T, nu, Cov_Mat] = ...
                    fun_MH_Candidate_Rob(accept_temp, k, W_n_LT, W_n_thresh, W_n_Cov_Mat, CrossFeed_Mat_Temp_cand, CrossFeed_Mat_Temp, Threshold_CF_cand, Threshold_CF, Struct_CrossFeed_Mat, p, n, liste_nrj_tot, liste_nrj_CF,...
                    acceptance_ratio_tot, acceptance_ratio_CF, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, theta_T, nu, Cov_Mat, Cov_Lt_eye, max_val(:,1), nb_rep, LT_CF_Death, LT_CF_Death, Weight_species, weight_day_seven, Temp);
        
                CF_Trace_Plot(:, :, n) = CrossFeed_Mat_Temp;  
        
                CrossFeed_Mat_Temp_cand(1:1+size(CrossFeed_Mat_Temp_cand,1):end) = 0;
                CrossFeed_Mat_Temp_cand(CrossFeed_Mat_Temp_cand < 1e-06) = 0;

                CrossFeed_Mat_Temp = CrossFeed_Mat_Temp.*(mu_max_vect_temp);
                CrossFeed_Mat_Temp_cand = CrossFeed_Mat_Temp_cand.*(mu_max_vect_temp);
                Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp_cand./mu_max_vect_temp;
            end
            Temp = Temp*0.99;

            Struct_CrossFeed_Mat{n + 1} = CrossFeed_Mat_Temp;
            n = n + 1;
            p = p + 1;
        
    
            if mod(p, num_steps) == 0
                disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
            end          
        end
        sum(accept_temp)/num_steps
        Cov_Mat_CF = Cov_Mat;
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
                    Threshold_CF, Threshold_death_cand, Threshold_Pred, Death_Mat_Temp_cand, death_rate,...
                    Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
                z_temp = deval(sol, Time_step);
                X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
                P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
                X = X + P;
                nu = k^(-3/4); 
                StackPlotTot = X./sum(X);
        
                [accept_temp, W_n_thresh_death, W_n_Cov_Mat_death, Death_Mat_Temp_cand, Death_Mat_Temp, Threshold_death_cand, Threshold_death, Struct_Death_Mat, liste_nrj_tot, liste_nrj_pred, acceptance_ratio_tot, acceptance_ratio_pred, p, r, theta_T_death, ~, nu, Cov_Mat] = ...
                    fun_MH_Candidate_Rob_Death(accept_temp, k, W_n_thresh_death, W_n_Cov_Mat_death, Death_Mat_Temp_cand, Death_Mat_Temp, Threshold_death_cand, Threshold_death, Struct_Death_Mat, p, r, liste_nrj_tot, liste_nrj_pred,...
                    acceptance_ratio_tot, acceptance_ratio_pred, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, theta_T_death, 0, nu, Cov_Mat, Cov_Lt_eye, max_val(:,2), nb_rep, LT_CF_Death, LT_CF_Death, Weight_species, weight_day_seven, Temp);
    
                Death_Trace_Plot(:,:, r) = Death_Mat_Temp;
                ind = eye(S,S);
                Death_Mat_Temp_cand(~ind) = 0; %Changement pour diagonale zeros
            end
            Temp = Temp*0.99;

            Struct_Death_Mat{r + 1} = Death_Mat_Temp;
            r = r + 1;
            p = p + 1;
    
            if mod(p, num_steps) == 0
                disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
            end
        end
        sum(accept_temp)/num_steps
        Cov_Mat_Death = Cov_Mat;
    
        %%%%%%%%%% Resource matrix loop %%%%%%%%
    
        Cov_Mat = 0.1*eye(S*nb_Res); %Cov_Mat_Res; %
        Cov_Mat_LT = 0.01*eye(S*nb_Res);
        W_n_LT_Res = normrnd(zeros(S*nb_Res,1), ones(S*nb_Res,1)); W_n_thresh_Res = normrnd(zeros(S,1), ones(S,1)); W_n_Cov_Mat_Res = normrnd(zeros(S*nb_Res,1), ones(S*nb_Res,1));
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

                [accept_temp, W_n_LT_Res, W_n_thresh_Res, W_n_Cov_Mat_Res, Resource_Matrix_cand, Resource_Matrix, T_cand, T, Struct_Res_Mat, liste_nrj_tot, liste_nrj_res, acceptance_ratio_tot, acceptance_ratio_res, p, l, ~, nu, Cov_Mat, Cov_Mat_LT, Lag_time_Cons_cand, Lag_time_Cons] = ...
                    fun_MH_Candidate_Rob(accept_temp, k, W_n_LT_Res, W_n_thresh_Res, W_n_Cov_Mat_Res, Resource_Matrix_cand, Resource_Matrix, zeros(S,1), zeros(S,1), Struct_Res_Mat, p, l, liste_nrj_tot, liste_nrj_res,...
                    acceptance_ratio_tot, acceptance_ratio_res, S, nb_Res, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, eye(S), nu, Cov_Mat, Cov_Mat_LT, max_val(:,3), nb_rep, Lag_time_Cons_cand, Lag_time_Cons, Weight_species, weight_day_seven, Temp);
            end
            Temp = Temp*0.99;

            Struct_Res_Mat{r + 1} = Resource_Matrix;
            l = l + 1;
            p = p + 1;


            if mod(p, num_steps) == 0
                disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
            end
        end
        sum(accept_temp)/(num_steps)
        Cov_Mat_Res = Cov_Mat;

        %%%%%%%%%% Predation matrix loop %%%%%%%%
        %If kept modify the threshold too

        Cov_Mat = 0.1*eye(S*S);
        Cov_Mat_LT = 0.01*eye(S*S);
        W_n_LT_Pred = normrnd(zeros(S*S,1), ones(S*S,1)); W_n_thresh_Pred = normrnd(zeros(S,1), ones(S,1)); W_n_Cov_Mat_Pred = normrnd(zeros(S*S,1), ones(S*S,1));
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

                [accept_temp, W_n_LT_Pred, W_n_thresh_Pred, W_n_Cov_Mat_Pred, Pred_Mat_Lyso_cand, Pred_Mat_Lyso, T_cand, T, Struct_Res_Mat, liste_nrj_tot, liste_nrj_res, acceptance_ratio_tot, acceptance_ratio_res, p, l, ~, nu, Cov_Mat] = ...
                    fun_MH_Candidate_Rob(accept_temp, k, W_n_LT_Pred, W_n_thresh_Pred, W_n_Cov_Mat_Pred, Pred_Mat_Lyso_cand, Pred_Mat_Lyso, zeros(S,1), zeros(S,1), Struct_Res_Mat, p, l, liste_nrj_tot, liste_nrj_res,...
                    acceptance_ratio_tot, acceptance_ratio_res, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, eye(S), nu, Cov_Mat, Cov_Mat_LT, max_val(:,3), nb_rep, LT_CF_Death, LT_CF_Death, Weight_species, weight_day_seven, Temp);

                Pred_Mat_Lyso_cand = Binary_Pred.*Pred_Mat_Lyso_cand; %multiplication by the binary predation matrix to get a 0-1 matrix. 1 at imposed prey 0 otherwise. We don't seek the preys, we assume there are known. To discuss but otherwise I'm afraid there would have too many possibilites (already true)
            end
            Temp = Temp*0.99;

            Struct_Pred_Mat_Lyso{r + 1} = Pred_Mat_Lyso;
            h = h + 1;
            p = p + 1;


            if mod(p, num_steps) == 0
                disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
            end
        end
        sum(accept_temp)/(num_steps)
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
    
    % [Fin_mat_1, Sigma_mat_1, h_mat_1, p_mat_1, Slope_mat_1, Intercept_mat_1] = Val_Trace(CF_Trace_Plot, 100);
    % CrossFeed_Mat_Temp(h_mat_1 == 0) = Fin_mat_1(h_mat_1 == 0);
    % [Fin_mat, Sigma_mat, h_mat, p_mat, Slope_mat, Intercept_mat] = Val_Trace(Res_Trace_Plot, 100);
    % Resource_Matrix(h_mat == 0) = Fin_mat(h_mat == 0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Data to test on half of the remaining data
    Time_step = [0 1 3 7 10 21]*24; %Measured time step in hours
    tspan = [0, max(Time_step)]; %Time interval in hours
    Data_Evol_temp = table2array(Data_Evol_test(:, 2:end)); %table2array(Data_Evol(:, 2:end));
    % Data_Evol_temp = Data_Evol_temp(:,ismember(mod(1:length(Data_Evol_temp(1,:)), 4), [1, 2, 3]));%20% of the data to test
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
    
    %CrossFeed_Mat_Temp_rates = CrossFeed_Mat_Temp.*mu_max_vect_temp; Resource_Matrix_rates = Resource_Matrix.*mu_max_vect_temp;
            
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

    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).CrossFeed_Mat_Temp = CrossFeed_Mat_Temp;   
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Pred_Mat_Lyso = Pred_Mat_Lyso;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Threshold_Pred = Threshold_Pred;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Threshold_CF = Threshold_CF;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Threshold_death = Threshold_death;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Death_Mat_Temp = Death_Mat_Temp;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Resource_Matrix = Resource_Matrix;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Lag_time_Pred = Lag_time_Pred;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Lag_time_Cons = Lag_time_Cons;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Mat_kappa_3 = Mat_kappa_3;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).kappa_mat = kappa_mat;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).R = R;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).CrossFeed_Mat_init = CrossFeed_Mat_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Threshold_CF_init = Threshold_CF_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Pred_Mat_Lyso_init = Pred_Mat_Lyso_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Threshold_death_init = Threshold_death_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Death_Mat_init = Death_Mat_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Resource_Matrix_init = Resource_Matrix_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Lag_time_Pred_init = Lag_time_Pred_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Lag_time_Cons_init = Lag_time_Cons_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Mat_kappa_3_init = Mat_kappa_3_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).kappa_mat_init = kappa_mat_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).R_init = R_init;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Struct_Res_Mat = Struct_Res_Mat;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Struct_Death_Mat = Struct_Death_Mat;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Struct_CrossFeed_Mat = Struct_CrossFeed_Mat;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Struct_Pred_Mat_Lyso = Struct_Pred_Mat_Lyso;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Cov_Mat_CF = Cov_Mat_CF;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Cov_Mat_Death = Cov_Mat_Death;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Cov_Mat_Res = Cov_Mat_Res;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).Cov_Mat_LT = Cov_Mat_LT;
    data_to_save_SynCom21_Soil_with_pred_Updated_Parfor(nb_iter_tot).death_rate = death_rate;
end

save(strcat(cd, '/Data/', 'data_to_save_SynCom21_Soil_with_pred_Updated_Parfor'), 'data_to_save_SynCom21_Soil_with_pred_Updated_Parfor')