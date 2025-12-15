%Simulations subsytems with fitted interspecific interactions
clear
close all


%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file = 'Death_CF_Model_New_Mu_max_v20';%'Death_CF_Model_V17';%

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','94:114', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');%readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');%
Data_Evol_test = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
S = height(Data_Evol);

%Load parameters
kappa_mat = load(strcat('Data/', Name_file, '_Kappa_mat.mat')); kappa_mat = kappa_mat.kappa_mat;
CrossFeed_Mat = load(strcat('Data/', Name_file, '_CrossFeed_Mat.mat')); CrossFeed_Mat = CrossFeed_Mat.CrossFeed_Mat_Temp;
Resource_Matrix = load(strcat('Data/', Name_file, '_Resource_Matrix.mat')); Resource_Matrix = Resource_Matrix.Resource_Matrix;
Pred_Mat = load(strcat('Data/', Name_file, '_Pred_Mat.mat')); Pred_Mat = Pred_Mat.Death_Mat_Temp;
Threshold = load(strcat('Data/', Name_file, '_Threshold.mat')); Threshold = Threshold.Threshold_CF;%Threshold.Threshold; For Philip data
Threshold_Pred = load(strcat('Data/', Name_file, '_Threshold_Pred.mat')); Threshold_Pred = Threshold_Pred.Threshold_death;
Lag_time_Cons = load(strcat('Data/', Name_file, '_Lag_time_Cons.mat')); Lag_time_Cons = Lag_time_Cons.Lag_time_Cons;
Lag_time_Pred = load(strcat('Data/', Name_file, '_Lag_time_Pred.mat')); Lag_time_Pred = Lag_time_Pred.Lag_time_Pred;
% death_rate = load(strcat('Data/', Name_file,'_death_rate.mat'), 'death_rate'); death_rate = death_rate.death_rate;
Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
R = load(strcat('Data/', Name_file, '_R_mat.mat')); R = R.R;
Name_file = 'Death_CF_Model_New_Mu_max_v19';%'Death_CF_Model_V17';
death_rate = load(strcat('Data/', Name_file,'_death_rate.mat'), 'death_rate'); death_rate = death_rate.death_rate;
name = string(table2array(Parameters_set(1:20,1)));
Time_step = [0 1 3 7 10 21]*24;%0:1:168;%Time step in hours

%Data to test on half of the remaining data
tspan = [0, max(Time_step)]; %Time interval in hours
Data_Evol_temp = table2array(Data_Evol_test(:, 2:end)); %table2array(Data_Evol(:, 2:end));
Data_Evol_temp = Data_Evol_temp(:,ismember(mod(1:length(Data_Evol_temp(1,:)), 4), [1, 2, 3]));%20% of the data to test
nb_time_step = length(Time_step);
nb_obs = length(Data_Evol_temp(1,:));
nb_rep = nb_obs/nb_time_step;
Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);

mean_y_0 = normrnd(mean(Measured_Abund(:,1,:), 3), std(Measured_Abund(:,1,:), 1, 3)); 
std_y_0 = std(Measured_Abund(:,1,:), 1, 3);


% Measured_Abund = table2array(Data_Evol(1:20, 2:7));
Threshold_Surviving = 1e-10;

Ratio_fin = zeros(1,20);
num_fig = 1;

diff_fin = zeros(1,S);

figure(num_fig); 
bar(mean(StackPlot_Meas,3)', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
num_fig = num_fig + 1;

nb_Res = length(R); %Number of resources (group of resources)
t_0 = 0; %Time 0
n_exp = 1; %No transfer in Philip experiment 
tspan = [0, max(Time_step)]; %[0, 22*24]; %Time interval in hours
S_sim = 1:S; %Species present into the subsystem
nb_iter = 1000;
diff_rel_biomass = zeros(1, nb_iter);
StackPlot_Sim_mean = zeros(S,6);
Abs_Abund_Sim_mean = zeros(S,6);
sum_biomass = sum(mean(Measured_Abund(:,1,:), 3));
for n = 1:nb_iter
    
    %Variation initial biomass
    mean_y_0 = lognrnd(log(mean(Measured_Abund(:,1,:), 3) + (5*std(Measured_Abund(:,1,:), 1, 3)).^2/2), (3e07*std(Measured_Abund(:,1,:), 1, 3))); 
    % mean_y_0 = normrnd(mean(Measured_Abund(:,1,:), 3), 2*ones(S,1)*mean(std(Measured_Abund(:,1,:), 1, 3))); 
    mean_y_0 = mean_y_0*sum_biomass/normrnd(sum(mean_y_0), 0.05*sum(mean_y_0)); %Variation in the total biomass
    
    kappa_mat_temp = kappa_mat(S_sim,:);
    CrossFeed_Mat_temp = CrossFeed_Mat(S_sim, S_sim);
    Mat_kappa_3_temp = Mat_kappa_3(S_sim, S_sim);
    Resource_Matrix_temp = Resource_Matrix(S_sim,:);
    Pred_Mat_temp = Pred_Mat(S_sim,S_sim);
    Threshold_temp = Threshold(S_sim);
    Threshold_Pred_temp = Threshold_Pred(S_sim);
    Lag_time_Cons_temp = Lag_time_Cons(S_sim,:);
    Lag_time_Pred_temp = Lag_time_Pred(S_sim,S_sim);
    name_temp = name(S_sim);
    Res_name = 1:nb_Res;
    yield_Pred_temp = 0;
    mean_y_0_temp = mean_y_0(S_sim);
    temp_bio = mean_y_0(S_sim);
    mean_y_0_temp(1) = 0*mean_y_0_temp(1);
    % temp_bio(n) = 100*mean(temp_bio);
    % mean_y_0_temp = temp_bio*sum(mean_y_0_temp)/sum(temp_bio);
    
    
    %Setting for Matlab ODE solver
    nb_tot_Species = S*3 + nb_Res;
    opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
    
    %Initialization of the colors
    
    colors = distinguishable_colors(60);
    
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0_temp;
    
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];

    y_0 = sum(mat_y_0(:,1:2),2);
    mat_y_0 = reshape(mat_y_0', 1, []);
    mat_y_0 = [mat_y_0 R];

    sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat_temp, CrossFeed_Mat_temp, Mat_kappa_3_temp, Resource_Matrix_temp, Threshold_temp, Threshold_Pred_temp, Pred_Mat_temp, death_rate, S, Lag_time_Cons_temp, Lag_time_Pred_temp, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); 
    z_temp = deval(sol, Time_step);
    X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
    W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
    R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass

    z_temp = z_temp(1:(end-nb_Res), end);
    z_temp = reshape(z_temp',3, S);
    z_temp = z_temp';
    z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
    StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
    [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
    % diff_biomass_temp = abs(Measured_Abund - X);
    StackPlotTot = X./sum(X);
    StackPlot_Sim_mean = StackPlot_Sim_mean + StackPlotTot;
    Abs_Abund_Sim_mean = Abs_Abund_Sim_mean + X;
    
    z_fin_sim = z(1:S,end); %Absolute stationary abundances
    nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
    z_fin_obs = mean(Measured_Abund(1:S, end, :), 3); %Absolute stationary abundances
    
    diff_fin(n) = sum(z_fin_obs - z_fin_sim);
    diff_rel_biomass_temp = StackPlot_Meas(:,2:end,:) - StackPlotTot(:,2:end);
    diff_rel_biomass(n) = sum(reshape(diff_rel_biomass_temp, 1, []));
end
StackPlot_Sim_mean = StackPlot_Sim_mean/nb_iter;
Abs_Abund_Sim_mean = Abs_Abund_Sim_mean/nb_iter;

figure(num_fig)
for j = 1:S
    plot(Time_step, Abs_Abund_Sim_mean(j,:), '-*', 'Color', colors(j,:));
    hold on
end
legend(name_temp, 'Orientation', 'vertical', 'Location', 'southeast')
num_fig = num_fig + 1;

hold on
for j = 1:S
    plot(Time_step, mean(Measured_Abund(j,:,:), 3), '--*', 'Color', colors(j,:));
    hold on
end

figure(num_fig);
bar(StackPlot_Sim_mean', 'stacked');
axis([0 11.5 0 1])
legend(name_temp, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;


figure(num_fig);
for j = 1:S
    scatter(Abs_Abund_Sim_mean(j, end), z_fin_obs(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors(j,:), 'MarkerFaceColor', colors(S+1-j,:));%col(S+1-j,:))      
    hold on
end
axis([0 3*10^(-4) 0 4*10^(-4)]);
reflin = refline(1,0);
axis square
reflin.Color = 'r';
xlabel('Experiment'); 
ylabel('Simulation');
legend(name);
title('Scatter absolute abundance');
num_fig = num_fig + 1;

[sorted_data, sorted_ind] = sort(diff_fin);
% sorted_name = name(sorted_ind);
figure(num_fig)
bar(sorted_ind, sorted_data)
num_fig = num_fig + 1;
figure(num_fig)
hist(diff_rel_biomass)
[h, p] = kstest((diff_rel_biomass - mean(diff_rel_biomass))/std(diff_rel_biomass));