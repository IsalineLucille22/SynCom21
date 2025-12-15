%Simulations subsytems with fitted interspecific interactions
clear
close all


%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file = 'Death_CF_Model_New_Mu_max_v28';%'Death_CF_Model_V17';%

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
death_rate = load(strcat('Data/', Name_file,'_death_rate.mat'), 'death_rate'); death_rate = death_rate.death_rate;
Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
R = load(strcat('Data/', Name_file, '_R_mat.mat')); R = R.R;
name = string(table2array(Parameters_set(1:20,1)));
Time_step = [0 1 3 7 10 21]*24;%0:1:168;%Time step in hours
nb_Res = length(R); %Number of resources (group of resources)

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

mean_y_0 = mean(Measured_Abund(:,1,:), 3);
std_y_0 = std(Measured_Abund(:,1,:), 1, 3);

Ratio_fin = zeros(1,20);
num_fig = 1;

diff_fin = zeros(1,S);
mat_y_0 = [mean_y_0 zeros(S,1) zeros(S,1)];
y_0 = sum(mat_y_0(:,1:2),2);
mat_y_0 = reshape(mat_y_0', 1, []);
mat_y_0 = [mat_y_0 R];

%Setting for Matlab ODE solver
nb_tot_Species = S*3 + nb_Res;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.

% sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat_temp, CrossFeed_Mat_temp, Mat_kappa_3_temp, Resource_Matrix_temp, Threshold_temp, Threshold_Pred_temp, Pred_Mat_temp, death_rate, S, Lag_time_Cons_temp, Lag_time_Pred_temp, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
sol = ode45(@(t, y) Emergent_Encoding_fun(t, y, sum(R), [Pred_Mat(1,1); Pred_Mat(2,2)], [kappa_mat(1,2); kappa_mat(2,2)]), tspan,  mean_y_0(1:2), opts_1); 
z_temp = deval(sol, Time_step);
% X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
% W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
% R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass

plot(Time_step, z_temp, 'o-')