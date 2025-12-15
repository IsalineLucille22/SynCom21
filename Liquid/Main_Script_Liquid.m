%Simulations with fitted interspecific interactions
clear
close all

addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')


%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file = 'data_to_save_liquid_3.mat';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');%%readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8, 'Range','48:69', 'Format','auto');%
% Data_Evol = readtable(strcat('Data/Liquid_Data/','abund_cfu_se.xlsx'), 'Sheet', 3, 'Range','26:47', 'Format','auto');%Data in liquid from Clara's experiments
Data_Evol = readtable(strcat('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/Liquid_Data/','abund_cfu_se.xlsx'), 'Sheet', 3, 'Range','26:47', 'Format','auto');%Data in liquid from Clara's experiments
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
Time_step = [0 12 24 48 168]; %Measured time step in hours
tspan = [0, max(Time_step)]; %Time interval in hours
S = height(Data_Evol);
Time_step_BP = 1:5:168;
    
%Load parameters
load(strcat('Data/', Name_file));
data_to_save = data_to_save_liquid_3;
% Fitted parameters
ind_sim = 6; %1 to 6 data sets possible
kappa_mat = data_to_save(ind_sim).kappa_mat;
CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
Resource_Matrix = data_to_save(ind_sim).Resource_Matrix;
Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; 
Threshold_CF = data_to_save(ind_sim).Threshold_CF;
Threshold_death = data_to_save(ind_sim).Threshold_death; %Change it
Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
R = data_to_save(ind_sim).R;
nb_Res = length(R); %Number of resources (group of resources)
Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
death_rate =  data_to_save(ind_sim).death_rate;
name = string(table2array(Parameters_set(1:21,1)));
Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
Resource_Matrix = Resource_Matrix.*Var_Resource_Matrix;
%From Lysobacter
Pred_Mat_Lyso = data_to_save(ind_sim).Pred_Mat_Lyso;
Threshold_Pred = data_to_save(ind_sim).Threshold_Pred;

% %Random
% Name_file = 'data_liquid_rand.mat'; %Name of the random data sets (Containing 10 random data sets)
% load(strcat('Data/', Name_file));
% data_to_save = data_liquid_rand;
% % Fitted parameters
% ind_sim = 6; %1 to 6 data sets possible
% kappa_mat = data_to_save(ind_sim).kappa_mat_init;
% CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_init;
% Resource_Matrix = data_to_save(ind_sim).Resource_Matrix_init;
% Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_init; 
% Threshold_CF = data_to_save(ind_sim).Threshold_CF_init;
% Threshold_death = data_to_save(ind_sim).Threshold_death_init;
% Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons_init; 
% Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred_init;
% R = data_to_save(ind_sim).R_init;
% nb_Res = length(R); %Number of resources (group of resources)
% Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
% death_rate =  data_to_save(ind_sim).death_rate;
% name = string(table2array(Parameters_set(1:21,1)));
% Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
% Resource_Matrix = Resource_Matrix.*Var_Resource_Matrix;
% %From Lysobacter
% Pred_Mat_Lyso = data_to_save(ind_sim).Pred_Mat_Lyso;
% Threshold_Pred = data_to_save(ind_sim).Threshold_Pred;
% 
% %No cross-feeding
% Name_file = 'data_liquid_no_CF.mat'; %Name of the random data sets (Containing 10 random data sets)
% load(strcat('Data/', Name_file));
% data_to_save = data_liquid_no_CF;
% % Fitted parameters
% ind_sim = 6; %1 to 6 data sets possible
% kappa_mat = data_to_save(ind_sim).kappa_mat;
% CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
% Resource_Matrix = data_to_save(ind_sim).Resource_Matrix;
% Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; 
% Threshold_CF = data_to_save(ind_sim).Threshold_CF;
% Threshold_death = data_to_save(ind_sim).Threshold_death; %Change it
% Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
% Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
% R = data_to_save(ind_sim).R;
% nb_Res = length(R); %Number of resources (group of resources)
% Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
% death_rate =  data_to_save(ind_sim).death_rate;
% name = string(table2array(Parameters_set(1:21,1)));
% Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
% Resource_Matrix = Resource_Matrix.*Var_Resource_Matrix;
% %From Lysobacter
% Pred_Mat_Lyso = data_to_save(ind_sim).Pred_Mat_Lyso;
% Threshold_Pred = data_to_save(ind_sim).Threshold_Pred;

Data_Evol_temp = table2array(Data_Evol(:, 2:end));
std_y_0 = 0;
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
end    
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
rand = 0.5;%unifrnd(0,1);
mean_y_0 = Measured_Abund(:, 1, 2) +  (1 - rand)*Measured_Abund(:, 1, 1);
mean_y_0 = mean_y_0 + normrnd(0, 0*0.05*mean(mean_y_0));
Measured_Abund = mean(Measured_Abund,3);


% Measured_Abund = table2array(Data_Evol(1:20, 2:7));
%Number of surviving species after 8 weeks
Threshold_Surviving = 1e-10;
nb_Surv_Obs = sum(Measured_Abund > Threshold_Surviving);

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
n_exp = 1; %No transfer in Philip experiment 
yield_Pred = 0.2;% 20% of yield for predation

%Setting for Matlab ODE solver
nb_tot_Species = S*3 + nb_Res;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.

%Initialization of the colors

colors = distinguishable_colors(S);

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

        sol = ode45(@(t, y) fun_CF_Death_Pred(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix, Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate, Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
        z_temp = deval(sol, Time_step);
        sum(z_temp(:,1))
        sum(z_temp(:,4))
        sum(z_temp(:,end))
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Byproduct growth biomass
        W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
        % X = X + P;
        %z_byproduct 
        z_temp_BP = deval(sol, Time_step_BP);
        X_BP = z_temp_BP(mod(1:S*3,3) == 1,:); %Species'biomass
        P_BP = z_temp_BP(mod(1:S*3,3) == 2,:); %Byproduct growth biomass
        W_BP = z_temp_BP(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp_BP = z_temp_BP(S*3 + 1: end,:); %Byproducts'biomass

        figure(num_fig)
        for j = 1:S
            plot(Time_step, X(j,:), '-o', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
            hold on
        end
        %legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 168 0 4e-04])
        num_fig = num_fig + 1;
        figure(num_fig)
        for j = 1:S
            plot(Time_step, Measured_Abund(j,:), '--*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
            hold on
        end
        %legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 168 0 4e-04])

        num_fig = num_fig + 1;
        z_temp = z_temp(1:(end-nb_Res), end);
        z_temp = reshape(z_temp',3, S);
        z_temp = z_temp';
        z = z_temp(:,1);%sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
        StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
        [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
        [Shann_obs, Simp_obs] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));
    end
    StackPlotTot = X./sum(X);
end

K_S_Waste = (CrossFeed_Mat + Mat_kappa_3)./kappa_mat(:,1); %Division by columns
K_S_Waste(K_S_Waste==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
for i = 2:(length(Time_step_BP) - 1)
    temp = zeros(S);
    for j = 1:S
        temp(j,:) = (P_BP(j,i + 1) - P_BP(j,i)).*(X_BP(j, i)*CrossFeed_Mat(j,:).*(W_BP(:,i)'./(W_BP(:,i)' + K_S_Waste(j,:))))/sum((X_BP(j, i)*CrossFeed_Mat(j,:).*(W_BP(:,i)'./(W_BP(:,i)' + K_S_Waste(j,:)))));
    end
    figure(num_fig)
    num_fig = num_fig + 1;
    heatmap(temp)
end


%Percentage contribution
for i = 1:length(Time_step) 
    temp_2 = zeros(S);  % Initialize the SxS contribution matrix
    for j = 1:S  % Loop over consumers
        BP_contrib = X(j, i)*CrossFeed_Mat(j,:).*(W(:,i)'./(W(:,i)' + K_S_Waste(j,:)));
        denom = sum(BP_contrib);
        if denom == 0
            denom = 1e-10;  % Small epsilon to avoid NaN
        end
        temp_2(j,:) = BP_contrib/denom;
    end
    figure(num_fig)
    num_fig = num_fig + 1;
    heatmap(temp_2)
end



StackPlotTot = StackPlotTot./nb_replicates;
z_fin_sim = z(1:S,end); %Absolute stationary abundances
z_fin_obs = Measured_Abund(1:S, end); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
sum(z_fin_sim)/sum(z_fin_obs)

h_vect = zeros(1,S);
p_vect = zeros(1,S);
for i = 1:S
    [h_vect(i),p_vect(i)] = ttest2(StackPlot_Meas(i,:)', StackPlotTot(i,:)');
end
name_diff = name(logical(h_vect));

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

[diff_vect, I, fact] = Sort_diff(StackPlotTot(:,end),StackPlot_Meas(:,end));
name_order = name(I);
h = figure(num_fig);
stem(diff_vect);
xtickangle(90)
set(gca,'xtick',1:21,'xticklabel',name_order)
ylabel('Sim - Obs')
num_fig = num_fig + 1;

figure(num_fig);
for j = 1:S
    scatter(StackPlot_Meas(j, n_exp+1), StackPlotTot(j,n_exp+1), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors(j,:), 'MarkerFaceColor', colors(S+1-j,:));     
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
    scatter(z_fin_obs(j), z_fin_sim(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors(j,:), 'MarkerFaceColor', colors(S+1-j,:));     
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

figure(num_fig);
plot(1:length(Time_step), Simp_obs, 'b--o')
hold on
plot(1:length(Time_step), Simp_sim, 'r--o')
num_fig = num_fig + 1;

X = X'; %For R sript comparison