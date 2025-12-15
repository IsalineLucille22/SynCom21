%Simulations with fitted interspecific interactions
clear
close all


%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file = 'data_to_save_Inter_v2';%'data_to_save_with_corrections';%'data_to_save_5.mat';%'data_to_save.mat';%'Fixed_Res_12_Groups_v1'; %'Death_CF_Model_New_Mu_max_v10';%'Death_CF_Model_V17';%'Rand_Parameters';%
Name_file_Saved = 'Fixed_Res_12_Groups_v2';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46','Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
S = height(Data_Evol);

%Load parameters
load(strcat('Data/', Name_file));
data_to_save = data_to_save_Inter_v2;%data_to_save_with_corrections;%data_to_save_5;%data_to_save_4;
% Fitted parameters
ind_sim = 1; %3, 5, (7), 9, (10), 11, 13, (15), (17), 18
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
Time_step = [0 1 3 7 10 15 21 22]*24; %[0 1 3 7 10 21]*24; 
Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
Resource_Matrix = Resource_Matrix.*Var_Resource_Matrix;
Time_step_BP = 1:50:500;


%%% Simulation with the mean of multiple datasets
% values = {data_to_save_Inter_v2.CrossFeed_Mat_Temp}; % Collect the field values into a cell array
% tot_mat = values{1, 1};
% for i = 2:length(values)
%     tot_mat = tot_mat + values{1, i};
% end
% CrossFeed_Mat =  tot_mat/length(values);


% %%% Random parameters simulations (one structure). Uncomment it to generate system with
% %%% random interactions.
% Name_file = 'data_random.mat'; %Name of the random data sets (Containing 10 random data sets)
% load(strcat('Data/', Name_file));
% data_to_save = data_random;
% ind_sim = 8; %from 1:10
% kappa_mat = data_to_save(ind_sim).kappa_mat_init;
% CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_init;
% Resource_Matrix = data_to_save(ind_sim).Resource_Matrix_init;
% Pred_Mat = data_to_save(ind_sim).Death_Mat_init; 
% Threshold = data_to_save(ind_sim).Threshold_CF_init;
% Threshold_Pred = data_to_save(ind_sim).Threshold_death_init;
% Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons_init; 
% Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred_init;
% R = data_to_save(ind_sim).R_init;
% nb_Res = length(R); %Number of resources (group of resources)
% Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
% death_rate =  data_to_save(ind_sim).death_rate;
% name = string(table2array(Parameters_set(1:20,1)));
% Time_step = [0 1 3 7 10 21]*24; 
% Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
% Resource_Matrix = Resource_Matrix.*Var_Resource_Matrix;
% 
%%% No cross-feeding parameters simulations (two structures, data_no_CF and data_no_CF_2 containing 10 sets). Uncomment it to generate
%%% system where cross-feeding interactions equal 0.
% Name_file = 'data_no_CF_2.mat'; %'data_no_CF.mat'; %Name of data sets with no cross-feeding (Containing 10 random data sets)
% load(strcat('Data/', Name_file));
% data_to_save = data_no_CF_2; %data_no_CF;
% % Fitted parameters
% ind_sim = 8; %from 1:10
% kappa_mat = data_to_save(ind_sim).kappa_mat;
% CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
% Resource_Matrix = data_to_save(ind_sim).Resource_Matrix;
% Pred_Mat = data_to_save(ind_sim).Death_Mat_Temp; 
% Threshold = data_to_save(ind_sim).Threshold_CF;
% Threshold_Pred = data_to_save(ind_sim).Threshold_death;
% Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
% Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
% R = data_to_save(ind_sim).R;
% nb_Res = length(R); %Number of resources (group of resources)
% Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
% death_rate =  data_to_save(ind_sim).death_rate;
% name = string(table2array(Parameters_set(1:20,1)));
% Time_step = [0 1 3 7 10 21]*24; 
% Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
% Resource_Matrix = Resource_Matrix.*Var_Resource_Matrix;

correction_table = [1,0.1,1,0.1,0.1,0.2,1,0.01,0.1,1,1,0.1,0.1,1,1,1,0.1,1,1,1]'; 
Data_Evol_temp = table2array(Data_Evol(:, 2:end));
std_y_0 = 0;
Data_Evol_temp = Data_Evol_temp(:, ismember(mod(1:length(Data_Evol_temp(1,:)),4), [1, 2, 3]));%20% of the data to train. Hold-out method.
% Data_Evol_temp(:, 1:3) = Data_Evol_temp(:, 1:3).*correction_table;
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


% Data_Evol_temp = table2array(Data_Evol(:, 2:end));
% std_y_0 = 0;
% Data_Evol_temp = Data_Evol_temp(:, ismember(mod(1:length(Data_Evol_temp(1,:)),4), [1, 2, 3]));%20% of the data to train. Hold-out method.
% % Data_Evol_temp(:, 1:4) = correction_table.*Data_Evol_temp(:, 1:4);%Only the inoculum is corrected
% nb_obs = length(Data_Evol_temp(1,:));
% nb_time_step = length(Time_step);
% nb_rep = nb_obs/nb_time_step;
% Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
% for i = 1:nb_rep
%     Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
% end
% StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
% rand = 0.5;%unifrnd(0,1);
% mean_y_0 = Measured_Abund(:, 1, 2) +  (1 - rand)*Measured_Abund(:, 1, 1);
% mean_y_0 = mean_y_0 + normrnd(0, 0*0.05*mean(mean_y_0));


% Measured_Abund = table2array(Data_Evol(1:20, 2:7));
%Number of surviving species after 8 weeks
Threshold_Surviving = 1e-10;
nb_Surv_Obs = sum(Measured_Abund > Threshold_Surviving);

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
n_exp = 1; %No transfer in Philip experiment 
tspan = [0, max(Time_step)*24]; %Time interval in hours
yield_Pred = 0;% 20% of yield for predation

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
            plot(Time_step, X(j,:), '-o', 'Color', colors{j});%, Time_step, R_temp, 'o');
            hold on
        end
        %legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 600 0 4e-04])
        num_fig = num_fig + 1;
        figure(num_fig)
        for j = 1:S
            plot(Time_step, Measured_Abund(j,:), '-*', 'Color', colors{j});%, Time_step, R_temp, 'o');
            hold on
        end
        legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 600 0 4e-04])

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

%Absolute contribution (estimated)
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
    title(strcat('Time point contribution', {' '}, num2str(Time_step_BP(i))))
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
    title(strcat('Time point contribution', {' '}, num2str(Time_step(i))))
end

StackPlotTot = StackPlotTot./nb_replicates;
z_fin_sim = z(1:S,end); %Absolute stationary abundances
z_fin_sim(15)
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
    scatter(StackPlot_Meas(j, n_exp+1), StackPlotTot(j,n_exp+1), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
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
    scatter(z_fin_obs(j), z_fin_sim(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
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

% %% Computation of the rates distribution for each species
% 
% [~, nb_iter] = size(data_to_save_Inter_v2);
% Mat_Struct = {data_to_save_Inter_v2.CrossFeed_Mat_Temp};
% [Sample_rate, LN_fit, nb_neg_val] = Rates_Distribution(12, Mat_Struct, num_fig, name);