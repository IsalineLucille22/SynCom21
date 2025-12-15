%Simulations subsytems with fitted interspecific interactions
clear
close all


%Save or Not
save_data = 0; %1 if save, 0 otherwise

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46', 'Format','auto');
Dry_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 5,'Range','29:49', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','94:114', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');%readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');%
Data_Evol_test = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
S = height(Data_Evol);
Dry_set = table2array(Dry_set(:, 4));


%Load parameters
yield_vect = table2array(Parameters_set(1:S,7)); %Change it according to the desired mu_max
Time_step = 0:1:504;%[0 1 3 7 10 21]*24;%0:1:168;%Time step in hours [0 1 3 7 10 21 31 61 92 122];%
Name_file = 'data_to_save_Inter_v2';%'data_to_save_5';%'data_to_save';%'data_to_save_Inter_v2';%
load(strcat('Data/', Name_file));
data_to_save = data_to_save_Inter_v2;%data_to_save_5;%data_to_save_4;
ind_sim = 3;
kappa_mat = data_to_save(ind_sim).kappa_mat;
CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
Resource_Matrix = data_to_save(ind_sim).Resource_Matrix; %repmat(yield_vect, 1, 12); %
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

%kappa_mat(:,2)./(kappa_mat(:,2) + kappa_mat(:,3))

% kappa_mat(:,3) = kappa_mat(:,2)./yield_vect - kappa_mat(:,2);

Data_Evol_temp = table2array(Data_Evol(:, 2:end));
std_y_0 = 0;
Data_Evol_temp = Data_Evol_temp(:,mod(1:length(Data_Evol_temp(1,:)), 4) == 2);%20% of the data to train. Hold-out method.
mean_y_0 = Data_Evol_temp(:,1);

%Data to test on half of the remaining data
tspan = [0, max(Time_step)]; %Time interval in hours
Data_Evol_temp = table2array(Data_Evol_test(:, 2:end)); %table2array(Data_Evol(:, 2:end));
Data_Evol_temp = Data_Evol_temp(:,ismember(mod(1:length(Data_Evol_temp(1,:)), 4), [1, 2, 3]));%20% of the data to test
nb_time_step = length(Time_step);
nb_obs = length(Data_Evol_temp(1,:));
nb_rep = nb_obs/nb_time_step;
% Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
% for i = 1:nb_rep
%     Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
% end
% StackPlot_Meas = Measured_Abund./sum(Measured_Abund);

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
tspan = [0, max(Time_step)]; %[0, 22*24]; %Time interval in hours
num_fig = 1;
max_biomass = zeros(1, length(yield_vect)); end_biomass = max_biomass;
for i = 1:length(yield_vect)
    S_sim = i; %Species present into the subsystem
    S = length(S_sim); %Number of species
    
    kappa_mat_temp = kappa_mat(S_sim,:);
    CrossFeed_Mat_temp = CrossFeed_Mat(S_sim, S_sim);
    Mat_kappa_3_temp = Mat_kappa_3(S_sim, S_sim);
    Resource_Matrix_temp = Resource_Matrix(S_sim,:);
    Pred_Mat_temp = 0*Pred_Mat(S_sim,S_sim);
    Threshold_temp = Threshold(S_sim);
    Threshold_Pred_temp = Threshold_Pred(S_sim);
    Lag_time_Cons_temp = Lag_time_Cons(S_sim,:);
    Lag_time_Pred_temp = Lag_time_Pred(S_sim,S_sim);
    name_temp = name(S_sim);
    mean_y_0_temp = mean_y_0(S_sim);
    death_rate_temp = 0*death_rate(S_sim);
    
    %Setting for Matlab ODE solver
    nb_tot_Species = S*3 + nb_Res;
    opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
    
    %Initialization of the colors
    
    colors = distinguishable_colors(60);
    
    nb_replicates = 1;
    
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0_temp;
    
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];
    
    
    y_0 = sum(mat_y_0(:,1:2),2);
    mat_y_0 = reshape(mat_y_0', 1, []);
    mat_y_0 = [mat_y_0 R];
    
    sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat_temp, CrossFeed_Mat_temp, Mat_kappa_3_temp, Resource_Matrix_temp, Threshold_temp, Threshold_Pred_temp, Pred_Mat_temp, death_rate_temp, S, Lag_time_Cons_temp, Lag_time_Pred_temp, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); 
    z_temp = deval(sol, Time_step);
    X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
    W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
    R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
    display(X)
    fig = figure; 
    figure(num_fig)
    for j = 1:S
        plot(Time_step, X(j,:), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
        hold on
    end
    axis([0 500 0 2e-03])
    yield_fin = (X(end)/sum(R))/yield_vect(i);
    legend_labels = strcat(name_temp, ", (ratio yield = ", string(yield_fin), ")");
    legend(legend_labels, 'Orientation', 'vertical', 'Location', 'northeast')
    num_fig = num_fig + 1;
    
    % figure(num_fig)
    % plot(Time_step, sum(R_temp), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
    % num_fig = num_fig + 1;
    % 
    % figure(num_fig)
    % for j = 1:S
    %     plot(Time_step, W(j,:), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
    %     hold on
    % end
    % num_fig = num_fig + 1;
    
    
    z_temp = z_temp(1:(end-nb_Res), end);
    z_temp = reshape(z_temp',3, S);
    z_temp = z_temp';
    z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
    nb_cell_per_mL = log10(X/(10*Dry_set(i)*1e-15));%nb cells per gram dry soil
    max_biomass(i) = max(nb_cell_per_mL);%max(X);
    end_biomass(i) = nb_cell_per_mL(end);%X(end);
    
    iFolderName = fullfile(cd, 'Figures'); % More robust path handling
    FigName = fullfile(iFolderName, strcat(Name_file, '_Test_Yield_', name_temp, '.pdf')); % Correct file format
    
    % Save the figure
    saveas(fig, FigName);
% close all
end

name_sorted = flip([2, 12, 6, 1, 13, 8, 18, 7, 16, 3, 5, 15, 17, 4, 19, 20, 10, 11, 9 14]);

figure(num_fig)
barh(max_biomass(name_sorted))
xlim([6 9])
yticks(1:length(name(name_sorted))); 
yticklabels(name(name_sorted));
num_fig = num_fig + 1;

figure(num_fig)
barh(end_biomass(name_sorted))
xlim([6 9])
yticks(1:length(name(name_sorted))); 
yticklabels(name(name_sorted));
num_fig = num_fig + 1;