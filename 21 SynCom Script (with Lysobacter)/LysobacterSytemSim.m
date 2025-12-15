%Script made by Adline Vouillamoz and Isaline Guex
clear;
close all;

%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file_Lyso = 'Rand_Lyso_2';

%%Loading data
cd ..
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','142:163', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto');
Data_Evol_test = Data_Evol;
Time_step = [0 1 3 7 10 21]*24; %Measured time step in hours
tspan = [0, max(Time_step)]; %Time interval in hours
S = height(Data_Evol);
name = string(table2array(Parameters_set(1:S,1)));

%%%Random parameters
yield_vect = table2array(Parameters_set(1:S,5))'; 
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8)); 
mu_max_vect = mu_max_dist(:,1);
Parameters_Senka_Lag_time = table2array(Parameters_Senka_Lag_time(:,7:8));
mean_param_LN = mu_max_dist(:,1);
var_param_LN = mu_max_dist(:,2);
mu_max_vect_temp = max(lognrnd(mu_max_vect, 0*mu_max_dist(:,2)), 0.001);
yield = max(normrnd(yield_vect, 0.0), 1e-05);
ind_Lyso = 6;

% % Load parameters
Name_file_load = 'Death_CF_Model_V9';%'Death_CF_Model'; %'Death_Model_LT'; %Name of the load data for initial setting
kappa_mat = load(strcat('Data/', Name_file_load, '_Kappa_mat.mat')); kappa_mat = kappa_mat.kappa_mat;
CrossFeed_Mat = load(strcat('Data/', Name_file_load, '_CrossFeed_Mat.mat')); CrossFeed_Mat = CrossFeed_Mat.CrossFeed_Mat_Temp;
Resource_Matrix = load(strcat('Data/', Name_file_load, '_Resource_Matrix.mat')); Resource_Matrix = Resource_Matrix.Resource_Matrix;
Death_Mat = load(strcat('Data/', Name_file_load, '_Pred_Mat.mat')); Death_Mat = Death_Mat.Death_Mat_Temp;%Pred_Mat_Temp
Threshold_CF = load(strcat('Data/', Name_file_load, '_Threshold.mat')); Threshold_CF = Threshold_CF.Threshold_CF;
Threshold_death = load(strcat('Data/', Name_file_load, '_Threshold_Pred.mat')); Threshold_death = Threshold_death.Threshold_death;%Threshold_pred.Threshold_pred;
Lag_time_Cons = load(strcat('Data/', Name_file_load, '_Lag_time_Cons.mat')); Lag_time_Cons = Lag_time_Cons.Lag_time_Cons;
Lag_time_Pred = load(strcat('Data/', Name_file_load, '_Lag_time_Pred.mat')); Lag_time_Pred = Lag_time_Pred.Lag_time_Pred;
R = load(strcat('Data/', Name_file_load, '_R_mat.mat')); R = R.R;
nb_Res = length(R);

% % %Random parameters
% kappa_mat_Lyso = Increased_Mat(kappa_mat, ind_Lyso, S, 4, 0);
% kappa_mat_Lyso(ind_Lyso,:) = [2.5e+05 mu_max_vect_temp(6) (mu_max_vect_temp(ind_Lyso)/yield(ind_Lyso) - mu_max_vect_temp(ind_Lyso)) 2.5e+05];
% Death_Mat_Lyso = Increased_Mat(Death_Mat, ind_Lyso, S, S, 1); 
% Death_Mat_Lyso(ind_Lyso, ind_Lyso) = lognrnd(log(0.02), mean(var_param_LN));
% CrossFeed_Mat_Lyso = Increased_Mat(CrossFeed_Mat, ind_Lyso, S, S, 1);
% CrossFeed_Mat_Lyso = CrossFeed_Mat_Lyso./kappa_mat_Lyso(:,2); %Normalized
% Lag_time_Cons_Lyso = Increased_Mat(Lag_time_Cons, ind_Lyso, S, nb_Res, 0);
% Lag_time_Cons_Lyso(ind_Lyso,:) = lognrnd(Parameters_Senka_Lag_time(ind_Lyso,1), 0*Parameters_Senka_Lag_time(ind_Lyso,2));
% x = rand(1, nb_Res);
% rand_indice = x > 0.3;
% Lag_time_Rand = normrnd(50, 10, 1, nb_Res);
% Lag_time_Cons_Lyso(ind_Lyso, rand_indice) = max(Lag_time_Cons_Lyso(rand_indice), Lag_time_Rand(rand_indice));
% Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
% Resource_Matrix_Lyso = Increased_Mat(Resource_Matrix, ind_Lyso, S, nb_Res, 0); 
% Resource_Matrix_Lyso = Resource_Matrix_Lyso.*Var_Resource_Matrix;
% Resource_Matrix_Lyso(ind_Lyso,:) = lognrnd(log(0.2), mean(var_param_LN), 1, nb_Res); 
% CrossFeed_Mat_Lyso(ind_Lyso, :) = lognrnd(log(0.2), mean(var_param_LN), 1, S);
% CrossFeed_Mat_Lyso(:, ind_Lyso) = lognrnd(log(0.2), mean(var_param_LN), S, 1);
% CrossFeed_Mat_Lyso(ind_Lyso, ind_Lyso) = 0;
% CrossFeed_Mat_Lyso = CrossFeed_Mat_Lyso.*kappa_mat_Lyso(:,2);
% Mat_kappa_3_Lyso = kappa_mat_Lyso(:,3).*CrossFeed_Mat_Lyso./kappa_mat_Lyso(:,2);
% Lag_time_Pred_Lyso = max(normrnd(150, 50, S, S), 0);
% Threshold_Pred = 1;
% Pred_Mat_Lyso = zeros(S,S);
% Prey_num = unique(randi(21, 10, 1));
% Prey_num(Prey_num == ind_Lyso) = [];
% Pred_Mat_Lyso(ind_Lyso, Prey_num) = lognrnd(log(0.1), mean(var_param_LN), 1, length(Prey_num));

% % % A priori
% kappa_mat_Lyso = Increased_Mat(kappa_mat, ind_Lyso, S, 4, 0);
% kappa_mat_Lyso(ind_Lyso,:) = [2.5e+05 mu_max_vect_temp(6) (mu_max_vect_temp(ind_Lyso)/yield(ind_Lyso) - mu_max_vect_temp(ind_Lyso)) 2.5e+05];
% Death_Mat_Lyso = Increased_Mat(Death_Mat, ind_Lyso, S, S, 1); 
% Death_Mat_Lyso(ind_Lyso, ind_Lyso) = mean([Death_Mat(5, 5), Death_Mat(13, 13)]);%lognrnd(log(0.02), mean(var_param_LN));%
% CrossFeed_Mat_Lyso = Increased_Mat(CrossFeed_Mat, ind_Lyso, S, S, 1);
% CrossFeed_Mat_Lyso = CrossFeed_Mat_Lyso./kappa_mat_Lyso(:,2); %Normalized
% Lag_time_Cons_Lyso = Increased_Mat(Lag_time_Cons, ind_Lyso, S, nb_Res, 0);
% Lag_time_Cons_Lyso(ind_Lyso,:) = mean(Lag_time_Cons_Lyso([5 13],:));%lognrnd(Parameters_Senka_Lag_time(ind_Lyso,1), 0*Parameters_Senka_Lag_time(ind_Lyso,2));
% % x = rand(1, nb_Res);
% % rand_indice = x > 0.3;
% % Lag_time_Rand = normrnd(50, 0, 1, nb_Res);
% % Lag_time_Cons_Lyso(ind_Lyso, rand_indice) = max(Lag_time_Cons_Lyso(rand_indice), Lag_time_Rand(rand_indice));
% Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
% Resource_Matrix_Lyso = Increased_Mat(Resource_Matrix, ind_Lyso, S, nb_Res, 0); 
% Resource_Matrix_Lyso = Resource_Matrix_Lyso.*Var_Resource_Matrix;
% Resource_Matrix_Lyso(ind_Lyso,:) = mean(Resource_Matrix_Lyso([5 13],:));%lognrnd(log(0.2), mean(var_param_LN), 1, nb_Res); %
% CrossFeed_Mat_Lyso(ind_Lyso, :) = mean(CrossFeed_Mat_Lyso([5 13], :));%lognrnd(log(0.2), mean(var_param_LN), 1, S);%
% CrossFeed_Mat_Lyso(:, ind_Lyso) = mean(CrossFeed_Mat_Lyso(:,[5 13]), 2);%lognrnd(log(0.2), mean(var_param_LN), S, 1);%
% CrossFeed_Mat_Lyso(ind_Lyso, ind_Lyso) = 0;
% CrossFeed_Mat_Lyso = CrossFeed_Mat_Lyso.*kappa_mat_Lyso(:,2);
% Mat_kappa_3_Lyso = kappa_mat_Lyso(:,3).*CrossFeed_Mat_Lyso./kappa_mat_Lyso(:,2);
% Lag_time_Pred_Lyso = max(normrnd(150, 0, S, S), 0);
% Threshold_Pred = 1;
% Pred_Mat_Lyso = zeros(S,S);
% Prey_num = [3 7 11 13 15 17 20];
% Pred_Mat_Lyso(ind_Lyso, Prey_num) = exp(mu_max_dist(ind_Lyso,1) + mu_max_dist(ind_Lyso,2).^2./2);%lognrnd(mu_max_dist(ind_Lyso,1), 0*mu_max_dist(ind_Lyso,2));%kappa_mat_Lyso(ind_Lyso,2)*ones(1, length(Prey_num));%Generate rates from fitted distribution

% %Fitted with Lyso
Name_file_load = 'Rand_Lyso_2';%'Lysobacter_Fitting_v3';%'Rand_Lyso';%
kappa_mat = load(strcat('Data/', Name_file_load, '_Kappa_mat.mat')); kappa_mat = kappa_mat.kappa_mat_Lyso;
CrossFeed_Mat = load(strcat('Data/', Name_file_load, '_CrossFeed_Mat.mat')); CrossFeed_Mat = CrossFeed_Mat.CrossFeed_Mat_Lyso;%CrossFeed_Mat.CrossFeed_Mat_Lyso_Temp;
Resource_Matrix = load(strcat('Data/', Name_file_load, '_Resource_Matrix.mat')); Resource_Matrix = Resource_Matrix.Resource_Matrix_Lyso;
Death_Mat = load(strcat('Data/', Name_file_load, '_Pred_Mat.mat')); Death_Mat = Death_Mat.Death_Mat_Lyso;%Death_Mat_Lyso_Temp;
Threshold_CF = load(strcat('Data/', Name_file_load, '_Threshold.mat')); Threshold_CF = Threshold_CF.Threshold_CF;
Threshold_death = load(strcat('Data/', Name_file_load, '_Threshold_death.mat')); Threshold_death = Threshold_death.Threshold_death;%'_Threshold_Pred.mat', '_Threshold_death.mat'
Threshold_Pred = load(strcat('Data/', Name_file_load, '_Threshold_Pred.mat')); Threshold_Pred = Threshold_Pred.Threshold_Pred;%'_Threshold_Pred_Lyso.mat'%
Lag_time_Cons = load(strcat('Data/', Name_file_load, '_Lag_time_Cons.mat')); Lag_time_Cons = Lag_time_Cons.Lag_time_Cons_Lyso;
Lag_time_Pred = load(strcat('Data/', Name_file_load, '_Lag_time_Pred.mat')); Lag_time_Pred = Lag_time_Pred.Lag_time_Pred_Lyso;
Pred_Mat_Lyso = load(strcat('Data/', Name_file_load, '_Pred_Mat_Lyso.mat')); Pred_Mat_Lyso = Pred_Mat_Lyso.Pred_Mat_Lyso;
R = load(strcat('Data/', Name_file_load, '_R_mat.mat')); R = R.R;
nb_Res = length(R);
kappa_mat_Lyso = kappa_mat;
Death_Mat_Lyso = Death_Mat; 
CrossFeed_Mat_Lyso = CrossFeed_Mat;
Lag_time_Cons_Lyso = Lag_time_Cons;
Resource_Matrix_Lyso = Resource_Matrix; 
Mat_kappa_3_Lyso = kappa_mat_Lyso(:,3).*CrossFeed_Mat_Lyso./kappa_mat_Lyso(:,2);
Lag_time_Pred_Lyso = Lag_time_Pred;

Data_Evol_temp = table2array(Data_Evol(:, 2:end));
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
rand = unifrnd(0, 1, 1, 3);
rand_1 = rand(1)/sum(rand); rand_2 = rand(2)/sum(rand); rand_3 = rand(3)/sum(rand);
mean_y_0 = rand_1*Measured_Abund(:,1,4) + rand_2*Measured_Abund(:,1,1) + rand_3*Measured_Abund(:,1,3);
mean_y_0 = mean_y_0 + normrnd(0, 0*mean(mean_y_0));
std_y_0 = std(Measured_Abund(:,1,:), 1, 3);
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
StackPlot_Meas = mean(StackPlot_Meas, 3);
Measured_Abund = mean(Measured_Abund,3);

% Measured_Abund = table2array(Data_Evol(1:20, 2:7));
%Number of surviving species after 8 weeks
Threshold_Surviving = 1e-10;
nb_Surv_Obs = sum(Measured_Abund > Threshold_Surviving);

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
n_exp = 1; %No transfer in Philip experiment 

yield_Pred = 0.2;%of yield for predation

%Setting for Matlab ODE solver
nb_tot_Species = S*3 + nb_Res;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.

%Initialization of the colors
colors_ind = [4, 1, 10, 15, 12, 11, 3, 8, 6, 20, 18, 19, 2, 5, 21, 13, 9, 14, 7, 16, 17];
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
% mean_y_0(6) = 0;
for i = 1:nb_replicates
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0;
    
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];

    for k = 1:n_exp
        y_0 = sum(mat_y_0(:,1:2),2);
        mat_y_0 = reshape(mat_y_0', 1, []);
        mat_y_0 = [mat_y_0 R];
   
%         sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat_Lyso, CrossFeed_Mat_Lyso, Mat_kappa_3_Lyso, Resource_Matrix_Lyso, Threshold_CF, Threshold_death, Death_Mat_Lyso, yield_Pred, S, Lag_time_Cons_Lyso, Lag_time_Pred_Lyso, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
        sol = ode45(@(t, y) fun_CF_Death_Pred(t, y, kappa_mat_Lyso, CrossFeed_Mat_Lyso, Mat_kappa_3_Lyso, Resource_Matrix_Lyso, Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Lyso, Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons_Lyso, Lag_time_Pred_Lyso, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
        X = X + P;
        figure(num_fig)
        for j = 1:S
            plot(Time_step, X(j,:), '-o', 'Color', colors{j});%, Time_step, R_temp, 'o');
            hold on
        end
        %legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 600 0 3e-04])
        num_fig = num_fig + 1;
        figure(num_fig)
        for j = 1:S
            plot(Time_step, Measured_Abund(j,:), '-*', 'Color', colors{j});%, Time_step, R_temp, 'o');
            hold on
        end
        %legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 600 0 3e-04])

        num_fig = num_fig + 1;
        z_temp = z_temp(1:(end-nb_Res), end);
        z_temp = reshape(z_temp',3, S);
        z_temp = z_temp';
        z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
        StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
        [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
        [Shann_obs, Simp_obs] = Shannon_Simpson_Indices(S,StackPlot_Meas);
    end
    StackPlotTot = X./sum(X);
end

StackPlotTot = StackPlotTot./nb_replicates;
z_fin_sim = z(1:S,end); %Absolute stationary abundances
z_fin_obs = Measured_Abund(1:S, end); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
sum(z_fin_sim)/sum(z_fin_obs)

h_vect = zeros(1,S);
p_vect = zeros(1,S);
h_vect_1 = zeros(1,S);
p_vect_1 = zeros(1,S);
Cos_sim = zeros(1,S);
for i = 1:S
    [h_vect(i), p_vect(i)] = ttest2(StackPlot_Meas(i,:)', StackPlotTot(i,:)');%ttest2(Measured_Abund(i,:)', X(i,:)');%
    [p_vect_1(i), h_vect_1(i)] = ranksum(Measured_Abund(i,:)', X(i,:)');%ttest2(Measured_Abund(i,:)', X(i,:)');%
    Cos_sim(i) = Measured_Abund(i,:)*X(i,:)'/(norm(Measured_Abund(i,:))*norm(X(i,:))); % Cosine similarity
end
name_diff = name(logical(h_vect_1));
[p_vect_2, h_vect_2] = hotell2(Measured_Abund, X);%

h_vect_comp = zeros(1,nb_time_step);
p_vect_comp = zeros(1,nb_time_step);
for i = 1:nb_time_step
    [p_vect_comp(i), h_vect_comp(i)] = ranksum(Measured_Abund(:,i), X(:,i));%ttest2(Measured_Abund(i,:)', X(i,:)');%
end

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
    scatter(StackPlot_Meas(j, n_exp+1), StackPlotTot(j,n_exp+1), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});     
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
    scatter(z_fin_obs(j), z_fin_sim(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});   
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
StackPlotTot_old = StackPlotTot;
StackPlot_Meas_old = StackPlot_Meas;
StackPlotTot(6:20,:) = StackPlotTot(7:21,:);
StackPlotTot(21,:) = StackPlotTot_old(6,:);
name_old = name;
name(6:20) = name(7:21);
name(21) = name_old(6);
StackPlotTot(21,:) = StackPlotTot_old(6,:);
bar(StackPlotTot', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;

figure(num_fig);
StackPlot_Meas(6:20,:) = StackPlot_Meas(7:21,:);
StackPlot_Meas(21,:) = StackPlot_Meas_old(6,:);
bar(StackPlot_Meas', 'stacked');
axis([0 11.5 0 1])
%legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
num_fig = num_fig + 1;

if save_data == 1
    save(strcat(cd, '/Data/', Name_file_Lyso,'_CrossFeed_Mat.mat'), 'CrossFeed_Mat_Lyso');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Threshold.mat'), 'Threshold_CF');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Threshold_death.mat'), 'Threshold_death');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Threshold_Pred.mat'), 'Threshold_Pred');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Pred_Mat.mat'), 'Death_Mat_Lyso');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Resource_Matrix.mat'), 'Resource_Matrix_Lyso');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Lag_time_Pred.mat'), 'Lag_time_Pred_Lyso');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Lag_time_Cons.mat'), 'Lag_time_Cons_Lyso');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Mat_kappa_3.mat'), 'Mat_kappa_3_Lyso');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Kappa_mat.mat'), 'kappa_mat_Lyso');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_Pred_Mat_Lyso.mat'), 'Pred_Mat_Lyso');
    save(strcat(cd, '/Data/', Name_file_Lyso,'_R_mat.mat'), 'R');
end