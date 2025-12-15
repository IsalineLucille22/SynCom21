%Simulations with fitted interspecific interactions. Script to simulate the
%development of the community over an extended period.
clear
close all


%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file = 'data_to_save_Inter';%'data_to_save.mat';%'Fixed_Res_12_Groups_v1'; %'Death_CF_Model_New_Mu_max_v10';%'Death_CF_Model_V17';%'Rand_Parameters';%
Name_file_Saved = 'Fixed_Res_12_Groups_v2';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46','Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
S = height(Data_Evol);

%Load parameters
load(strcat('Data/', Name_file));
data_to_save = data_to_save_Inter;%data_to_save_4;%
ind_sim = 3; %3, 5, (7), 9, (10), 11, 13, (15), (17), 18
kappa_mat = data_to_save(ind_sim).kappa_mat;
CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
Resource_Matrix = data_to_save(ind_sim).Resource_Matrix;
Pred_Mat = data_to_save(ind_sim).Death_Mat_Temp; 
Threshold = data_to_save(ind_sim).Threshold_CF;
Threshold_Pred = data_to_save(ind_sim).Threshold_death;
Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
R = data_to_save(ind_sim).R;
nb_Res = length(R); %Number of resources (group of resources)
Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
death_rate =  data_to_save(ind_sim).death_rate;
name = string(table2array(Parameters_set(1:20,1)));
Time_step = [0 1 3 7 10 21]*24; %[0 1 3 7 15 22]*24;%Measured time step in hours
% Time_step = [0 1 3 7 10 15 21 22]*24; %Measured time step in hours
Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
Resource_Matrix = Resource_Matrix.*Var_Resource_Matrix;
% Threshold(14) = 2*Threshold(14);
% death_rate(14) = 5*death_rate(14);



Data_Evol_temp = table2array(Data_Evol(:, 2:end));
std_y_0 = 0;
Data_Evol_temp = Data_Evol_temp(:, ismember(mod(1:length(Data_Evol_temp(1,:)),4), [1, 2, 3]));%20% of the data to train. Hold-out method.
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
rand = unifrnd(0,1);
mean_y_0 = Measured_Abund(:, 1, 2) +  (1 - rand)*Measured_Abund(:, 1, 1);
mean_y_0 = mean_y_0 + normrnd(0, 0.05*mean(mean_y_0));
Measured_Abund = mean(Measured_Abund,3);

% Time_step = 0:1:168;%55568;%Measured time step in hours

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
n_exp = 1; %No transfer in Philip experiment 
tspan = [0, max(Time_step)]; %Time interval in hours

%Setting for Matlab ODE solver
nb_tot_Species = S*3 + nb_Res;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.

%Initialization of the colors

%Initialization of the colors
colors_ind = [4, 1, 10, 15, 12, 3, 8, 6, 20, 18, 19, 2, 5, 21, 13, 9, 14, 7, 16, 17];
colors_init = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#B6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#01665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC', '#4393C3'};%distinguishable_colors(S);
colors = {};
for i = 1:S
    colors{i} = colors_init{colors_ind(i)};
end

nb_replicates = 1;

num_fig = 1;
for i = 1:nb_replicates
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0;
    
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];

    for k = 1:n_exp
        y_0 = sum(mat_y_0(:,1:2),2);
        mat_y_0 = reshape(mat_y_0', 1, []);
        mat_y_0 = [mat_y_0 R];

        sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix, Threshold, Threshold_Pred, Pred_Mat, death_rate, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        sum(z_temp(:,1))
        sum(z_temp(:,4))
        sum(z_temp(:,end))
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
        X = X + P;
        figure(num_fig)
        for j = 1:S
            plot(Time_step, X(j,:), '-', 'Color', colors{j});%, Time_step, R_temp, 'o');
            hold on
        end
        legend(name, 'Orientation', 'vertical', 'Location', 'southeast')

        num_fig = num_fig + 1;
        z_temp = z_temp(1:(end-nb_Res), end);
        z_temp = reshape(z_temp',3, S);
        z_temp = z_temp';
        z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
        StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
        [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
    end
    StackPlotTot = X./sum(X);
end

StackPlotTot = StackPlotTot./nb_replicates;
z_fin_sim = z(1:S,end); %Absolute stationary abundances

figure(num_fig);
bar(StackPlotTot', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;

figure(num_fig);
plot(1:length(Time_step), Simp_sim, 'r--o')
num_fig = num_fig + 1;