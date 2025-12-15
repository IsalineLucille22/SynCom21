%Script made by Adline Vouillamoz and Isaline Guex
clear;
close all;

%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file = 'Death_Model_Loop_2';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','94:114', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');
Data_Evol_test = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto');
S = height(Data_Evol);
name = string(table2array(Parameters_set(1:20,1)));

Set_CrossFeed_Mat = load(strcat('Data/', Name_file,'_Set_CrossFeed_Mat.mat')); Set_CrossFeed_Mat = Set_CrossFeed_Mat.Set_CrossFeed_Mat;
Set_T_CF = load(strcat('Data/', Name_file,'_Set_T_CF.mat')); Set_T_CF = Set_T_CF.Set_T_CF;
Threshold_death = load(strcat('Data/', Name_file,'_Threshold_Pred.mat')); Threshold_death = Threshold_death.Threshold_death;
Death_Mat_Temp = load(strcat('Data/', Name_file,'_Pred_Mat.mat')); Death_Mat_Temp = Death_Mat_Temp.Death_Mat_Temp;
Resource_Matrix = load(strcat('Data/', Name_file,'_Resource_Matrix.mat')); Resource_Matrix = Resource_Matrix.Resource_Matrix;
Lag_time_Pred = load(strcat('Data/', Name_file,'_Lag_time_Pred.mat')); Lag_time_Pred = Lag_time_Pred.Lag_time_Pred;
Lag_time_Cons = load(strcat('Data/', Name_file,'_Lag_time_Cons.mat')); Lag_time_Cons = Lag_time_Cons.Lag_time_Cons;
kappa_mat = load(strcat('Data/', Name_file,'_Kappa_mat.mat'), 'kappa_mat'); kappa_mat = kappa_mat.kappa_mat;
R = load(strcat('Data/', Name_file,'_R_mat.mat')); R = R.R;
yield_Pred = 0.2;
nb_Res = length(R);
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8)); %Normal Parameters_Senka_mu_max(:,2:3). Log-normal Parameters_Senka_mu_max(:,7:8)
mu_max_vect = mu_max_dist(:,1);
mu_max_vect_temp = exp(mu_max_vect + mu_max_dist(:,2).^2./2);

CrossFeed_Mat_Temp = Average_Mat(Set_CrossFeed_Mat); Threshold_CF = Average_Mat(Set_T_CF);
[test_dist, test_dist_rand, diff_obs, diff_rand, p_obs, p_rand] = Comparison_sample(Set_CrossFeed_Mat);
Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp./mu_max_vect_temp;

%Setting for Matlab ODE solver
%nb_tot_Species = S*3 + nb_Res;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
%Initialization of the colors
colors = distinguishable_colors(S);

%Data to test on half of the remaining data
Time_step = [0 1 3 7 10 21]*24; %Measured time step in hours
tspan = [0, max(Time_step)]; %Time interval in hours
Data_Evol_temp = table2array(Data_Evol_test(:, 2:end));
nb_time_step = length(Time_step);
nb_obs = length(Data_Evol_temp(1,:));
nb_rep = nb_obs/nb_time_step;
Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
gamma_Stack = std(StackPlot_Meas, 1, 3);
Threshold_Surviving = 1e-10;
nb_Surv_Obs = sum(mean(Measured_Abund,3) > Threshold_Surviving);
mean_y_0 = mean(Measured_Abund(:,1,:), 3);%table2array(Data_Evol(1:20, 2));
mat_y_0 = [mean_y_0 zeros(S,1) zeros(S,1)];
mat_y_0 = reshape(mat_y_0', 1, []);
mat_y_0 = [mat_y_0 R];

sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix, Threshold_CF, Threshold_death, Death_Mat_Temp, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); %Multiple resource groups
z_temp = deval(sol, Time_step);
X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass

CrossFeed_Mat_Temp_rates = CrossFeed_Mat_Temp.*mu_max_vect_temp; Resource_Matrix_rates = Resource_Matrix.*mu_max_vect_temp;
        
%%% Figure generation %%%
num_fig = 1;
ratio_fin_abund = sum(X(:,end))/sum(Measured_Abund(:,end));

figure(num_fig) %figure 1 : Simulated absolute abundances
for j = 1:S
    plot(Time_step, X(j,:), '-o', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
    hold on
end
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')

num_fig = num_fig + 1;
figure(num_fig) %figure 2 : Observed absolute abundances
for j = 1:S
    plot(Time_step, mean(Measured_Abund(j,:,:), 3), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
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
[Shann_obs, Simp_obs] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));
rand_prop = max(normrnd(1, 0.3, S, 1), 0);
mat_y_0 = (z_temp(:,1:2)/10).*rand_prop; %Take the 10% of the system
mat_y_0 = [mat_y_0 zeros(S,1)];

StackPlotTot = X./sum(X);

z_fin_sim = z(1:S,end); %Absolute stationary abundances
z_fin_obs = mean(Measured_Abund(1:S, end, :), 3); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
sum(z_fin_sim)/sum(z_fin_obs)

figure(num_fig); %figure 3 abondances relatives simulations
bar(StackPlotTot', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;

figure(num_fig); %figure 4 Relative observed abundances
bar(mean(StackPlot_Meas,3)', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
num_fig = num_fig + 1;

[diff_vect, I, ~] = Sort_diff(StackPlotTot(:,end),mean(StackPlot_Meas(:,end,:),3));
name_order = name(I);
figure(num_fig);  %figure 5 Difference between observed and simulated resluts at the end of the experiment
stem(diff_vect);
xtickangle(90)
set(gca,'xtick',1:21,'xticklabel',name_order)
ylabel('Sim - Obs')

num_fig = num_fig + 1;

figure(num_fig); %figure 7 
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


figure(num_fig); %figure 7
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

figure(num_fig); %figure 8
plot(1:length(Time_step), Simp_obs, 'b--o')
hold on
plot(1:length(Time_step), Simp_sim, 'r--o')
num_fig = num_fig + 1;

figure(num_fig); %figure 9
plot(1:length(Time_step), nb_Surv_Obs, 'b--o')
hold on
plot(1:length(Time_step), nb_Surv_Sim, 'r--o')
num_fig = num_fig + 1;