%Compare observed individual growth in soil of Syncom21 members with observed growth in Syncom20 mixture, and with modeled SynCom20 growth without cross-feeding.

%clearvars -except z 
clear
close all

%data path for individual growth
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')
name_Exp = 'Senka_21_CF_Liquid';

%% plot individual biomass in soil %Doesn't make sense here because not the same environment


%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Data_Evol = readtable(strcat('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/Liquid_Data/','abund_cfu_se.xlsx'), 'Sheet', 3, 'Range','26:47', 'Format','auto');%Data in liquid from Clara's experiments

biomass = readtable('MergedData_copy.xlsx','sheet','Growth kinetic params','Range','U3:X24','VariableNamingRule','preserve');
strains = readtable('MergedData_copy.xlsx','sheet','Growth kinetic params','Range','C3:C24','VariableNamingRule','preserve');

combined = [strains, biomass]; %Individual biomass in soil
combined = sortrows(combined,'Name','ascend');
combined = [combined, table(mean(table2array(combined(:,2:5)),2))];
combined = renamevars(combined,{'Var1'},{'mean'});
S = size(strains); S = S(1);


%actual absolute abundance data, has quadriplicate values for each time point

name = string(table2array(Parameters_set(1:S,1)));

% get single cell weight list

single_cell_mass = readtable(strcat('Data/','MergedData_copy.xlsx'), 'Sheet', 'Growth kinetic params','Range','C3:D24','Format','auto');
single_cell_mass = sortrows(single_cell_mass,'Name','ascend');

%sort alphabetically to strains
[name, idx2]= sortrows(name);

Data_Evol = Data_Evol(idx2,:);

%Automatized that
max_obs_biomass = mean(table2array(Data_Evol(:,18:S)),2); %timepoint 7

max_obs_biomass = [table(name), array2table(max_obs_biomass), Data_Evol(:,18:S)];%[table(name),array2table(max_obs_biomass), Data_Evol(:,22:25)];


%make a detour to have the strain names from top to bottom in alphabetical order

inv_names = sortrows(combined,'Name','descend');
xvalues = categorical(inv_names.Name);
X = reordercats(xvalues, inv_names.Name);
yvalues = max_obs_biomass.max_obs_biomass./(single_cell_mass.assumedCellDwInFgC*1e-15);

%make sure to flip the yvalues to have the inverted order

figH = figure;

subplot('Position',[0.3 0.1 0.2 0.62]);

barh(X,flip(log10(yvalues))) %observed biomass in Syncom21
hold on
barh(X,flip(log10(combined.mean)), 'FaceAlpha', 0.5)
%Scatter for each replicate
for i = 2:5
    scatter(flip(log10(table2array(combined(:,i)))),1:S,50,'.','r')
end
%Automatized this
for i = 3:6
    scatter(flip(log10(table2array(max_obs_biomass(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15)))),1:21,50,'.','b')
end

xlim([3 9]);
grid on
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'comparison-individual-SynCom21-biomass-in-soil.pdf');
saveas(figH,FigName,'pdf');

% ttest differences

ind_syncom_test = mattest(table2array(combined(:,2:5)),table2array(max_obs_biomass(:,3:6)./(single_cell_mass.assumedCellDwInFgC*1e-15)));
FDR = mafdr(ind_syncom_test, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing
ind_syncom_test=[combined.Name table(ind_syncom_test) table(FDR)];
writetable(ind_syncom_test, 'Data/individual-syncom-mattest.csv');

% take the difference of the means as fraction of the SynCom20 biomass

diff_vector = max_obs_biomass.max_obs_biomass-(combined.mean.*single_cell_mass.assumedCellDwInFgC*1e-15);
diff_fraction = diff_vector./max_obs_biomass.max_obs_biomass;

close all

figH = figure;
subplot('Position',[0.1 0.3 0.4 0.3]);
stem(diff_vector);
xtickangle(90)
set(gca,'xtick',1:S,'xticklabel',combined.Name)
ylabel('biomass difference');
title('Syncom - individual');
grid on
ylim([-4e-4 4e-4]);

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'syncom-individual-biomass-diff.pdf');
saveas(figH,FigName,'pdf');
%% comparison of Syncom actual growth to Syncom prediction with/without cross-feeding
% run script 'Plot_syncom_from_multi_simulation_data_no_CF'
% variable 'z' has the predicted species' biomass wiht/without cross-feeding 

load(strcat('Data/', name_Exp, 'z.mat'))
CF_biomass = z; %Reorder z according to the names
CF_biomass(CF_biomass < 0) = 0;
nb_data_set = size(CF_biomass, 2) - 1;
%the variable 'name' corresponds to the strain order on no_CF_biomass

CF_biomass = [table(name) array2table(CF_biomass)];
CF_biomass = sortrows(CF_biomass,'name','ascend');

no_CF_pop_counts=[table(name) array2table(max(z./(single_cell_mass.assumedCellDwInFgC*1e-15), 1))]; %Smaller than the dry biomass per gram per mL correct it
no_CF_pop_counts=sortrows(no_CF_pop_counts,'name','ascend');

% plot on top of each other

X = reordercats(xvalues,inv_names.Name);

close all

figH = figure;

subplot('Position',[0.3 0.1 0.2 0.62]);

barh(X,flip(log10(yvalues))) %observed biomass in Syncom20
hold on
barh(X,flip(log10(mean(table2array(no_CF_pop_counts(:, 2:(nb_data_set + 1))),2))), 'FaceAlpha', 0.5)
for i = 3:6
    scatter(flip(log10(table2array(max(max_obs_biomass(:,i)./(single_cell_mass.assumedCellDwInFgC*1e-15), 1)))),1:S,50,'.','b')
end
for i = 2:(nb_data_set + 1)
    scatter(flip(log10(table2array(no_CF_pop_counts(:,i)))),1:S,50,'.','r')
end

xlim([3 9]);
grid on
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'comparison-SynCom20-model-CF-biomass-in-soil.pdf');
saveas(figH,FigName,'pdf');

% ttest differences

syncom_CF_test = mattest(table2array(CF_biomass(:, 2:(nb_data_set + 1))),table2array(max_obs_biomass(:,3:6)));

FDR = mafdr(syncom_CF_test, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing

syncom_CF_test=[combined.Name table(syncom_CF_test) table(FDR)];

writetable(syncom_CF_test, 'Data/syncom-model-CF-mattest.csv');

% take the difference of the means as fraction of the SynCom20 biomass
diff_vector = max_obs_biomass.max_obs_biomass - mean(table2array(CF_biomass(:,2:(nb_data_set + 1))),2);

diff_fraction = diff_vector./max_obs_biomass.max_obs_biomass;

close all

figH = figure;

subplot('Position',[0.1 0.3 0.4 0.3]);

stem(diff_vector);
xtickangle(90)
set(gca,'xtick',1:S,'xticklabel',combined.Name)
ylabel('biomass difference');
title('Syncom - model-CF');
grid on
ylim([-4e-4 4e-4]);
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'syncom-model_CF-biomass-diff.pdf');
saveas(figH,FigName,'pdf');

%% Comparison individual Observed vs simulated. We simulate the biomass after one week and not 21 days. %Doesn't make sense here because not the same environment
close all 

num_fig = 1;
figH = figure;

subplot('Position',[0.3 0.1 0.2 0.62]);

barh(X,flip(log10(combined.mean)))
%Scatter for each replicate
hold on
for i=2:5
    scatter(flip(log10(table2array(combined(:,i)))),1:21,50,'.','r')
end

xlim([6 9]);
grid on
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'individual-SynCom21-observed-biomass-in-soil.pdf');
saveas(figH,FigName,'pdf');

%Simulations soil mono-cultures

Name_file =  'data_to_save_SynCom21_Soil_only_Lyso_Parfor_7';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_New';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Random';%'data_to_save_Inter_v2';%'data_no_CF.mat';
ind_sim = 1; %Chose the parameters set used


%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');
Dry_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 5,'Range','29:50', 'Format','auto');

%actual absolute abundance data, has quadriplicate values for each time point

%Data_Evol = readtable(strcat('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/Liquid_Data/','abund_cfu_IG.xlsx'), 'Sheet', 2, 'Range','46:67', 'Format','auto');%1:22 without correction for 0. Data in soil extract from Phil and Clara's experiments
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto'); %Senka's data
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
S = height(Data_Evol);
Time_step = [0 1 3 7 10 21]*24;%[0 12 48 96 7*24 21*24];%Measured time step in hours
yield_vect = table2array(Parameters_set(1:S,7)); %Change it according to the desired mu_max
Dry_set = table2array(Dry_set(:, 8));

Data_Evol_temp = table2array(Data_Evol(:, 2:end));
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
mean_y_0 = mean(Data_Evol_temp(:,1:nb_rep), 2); 

%Load parameters
load(strcat('Data/', Name_file));
data_to_save = data_to_save_SynCom21_Soil_only_Lyso_Parfor_7;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_Random;%data_to_save_Inter_v2;%data_no_CF;


%Initialization of the colors
colors_ind = [4, 1, 10, 15, 12, 3, 8, 6, 20, 18, 19, 2, 5, 21, 13, 9, 14, 7, 16, 17, 11];
colors_init = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#B6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#01665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC', '#4393C3', '#B35806'};%distinguishable_colors(S);
colors = {};
for i = 1:S
    colors{i} = colors_init{colors_ind(i)};
end

%Setting for Matlab ODE solver
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
nb_data_set = size(data_to_save, 2);

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
tspan = [0, max(Time_step)]; %[0, 22*24]; %Time interval in hours
num_fig = num_fig + 1;
max_biomass = zeros(length(yield_vect), nb_data_set); end_biomass = max_biomass;
yield_Pred = 0.2;
Time_step = [0 1 3 7]*24;%One-week time step
for j = 1:nb_data_set
    ind_sim = j;
    %Initialization of the interactions 
    
    kappa_mat = data_to_save(ind_sim).kappa_mat;
    CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
    Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; 
    Threshold_CF = data_to_save(ind_sim).Threshold_CF;
    Threshold_death = data_to_save(ind_sim).Threshold_death; 
    Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
    Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
    death_rate =  0*data_to_save(ind_sim).death_rate;
    name = string(table2array(Parameters_set(1:21,1)));
    Pred_Mat_Lyso = data_to_save(ind_sim).Pred_Mat_Lyso;
    Threshold_Pred = data_to_save(ind_sim).Threshold_Pred;
    R = data_to_save(ind_sim).R; 
    nb_Res = length(R); %Number of resources (group of resources)
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix; %Addition of a line for Lysobacter
    
    %Modifications upper part
    for i = 1:length(yield_vect)
        S_sim = i; %Species present into the subsystem
        S = length(S_sim); %Number of species
       
        kappa_mat_temp = kappa_mat(S_sim,:);
        CrossFeed_Mat_temp = CrossFeed_Mat(S_sim, S_sim);
        Mat_kappa_3_temp = Mat_kappa_3(S_sim, S_sim);
        Resource_Matrix_temp = Resource_Matrix(S_sim,:);
        Pred_Mat_Lyso_temp = 0*Pred_Mat_Lyso(S_sim,S_sim);
        Death_Mat_Temp_temp = Death_Mat_Temp(S_sim,S_sim);
        Threshold_CF_temp = Threshold_CF(S_sim);
        Threshold_Pred_temp = Threshold_Pred(S_sim);
        Threshold_death_temp = Threshold_death(S_sim);
        Lag_time_Cons_temp = Lag_time_Cons(S_sim,:);
        Lag_time_Pred_temp = Lag_time_Pred(S_sim,S_sim);
        name_temp = name(S_sim);
        mean_y_0_temp = mean_y_0(S_sim);
        death_rate_temp = death_rate(S_sim);
        
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
        
        sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat_temp, CrossFeed_Mat_temp, Mat_kappa_3_temp, Resource_Matrix_temp,...
        Threshold_CF_temp, Threshold_death_temp, Threshold_Pred_temp, Death_Mat_Temp_temp, death_rate_temp,...
        Pred_Mat_Lyso_temp, yield_Pred, S, Lag_time_Cons_temp, Lag_time_Pred_temp, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
    
        % fig = figure; 
        % figure(num_fig)
        % for j = 1:S
        %     plot(Time_step, X(j,:), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
        %     hold on
        % end
        % axis([0 500 0 2e-03])
        % yield_fin = (X(end)/sum(R))/yield_vect(i);
        % legend_labels = strcat(name_temp, ", (ratio yield = ", string(yield_fin), ")");
        % legend(legend_labels, 'Orientation', 'vertical', 'Location', 'northeast')
        % num_fig = num_fig + 1;
        
        z_temp = z_temp(1:(end-nb_Res), end);
        z_temp = reshape(z_temp',3, S);
        z_temp = z_temp';
        z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
        nb_cell_per_mL = log10(X/(10*Dry_set(i)*1e-15));%nb cells per gram dry soil
        max_biomass(i, j) = max(nb_cell_per_mL);%max(X);
        end_biomass(i, j) = nb_cell_per_mL(end);%X(end);
        
        % iFolderName = fullfile(cd, 'Figures'); % More robust path handling
        % FigName = fullfile(iFolderName, strcat(Name_file, '_Test_Yield_', name_temp, '.pdf')); % Correct file format
        
        % Save the figure
        % saveas(fig, FigName);
    % close all
    end
end

name_sorted = flip([2, 13, 7, 1, 14, 9, 19, 8, 17, 3, 6, 5, 16, 18, 4, 20, 21, 11, 12, 10, 15]);

% figure(num_fig)
% barh(max_biomass(name_sorted))
% xlim([6 9])
% yticks(1:length(name(name_sorted))); 
% yticklabels(name(name_sorted));
% num_fig = num_fig + 1;
% 
% figure(num_fig)
% barh(end_biomass(name_sorted))
% xlim([6 9])
% yticks(1:length(name(name_sorted))); 
% yticklabels(name(name_sorted));
% num_fig = num_fig + 1;

figG = figure;
subplot('Position',[0.3 0.1 0.2 0.62]);
mean_end_biomass = mean(end_biomass, 2);

barh(name(name_sorted), mean_end_biomass(name_sorted))
% Scatter for each replicate
hold on
for i = 1:nb_data_set
    scatter(end_biomass(name_sorted,i), 1:S_sim, 50, '.', 'r')
end
xlim([6 9]);
grid on
iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'individual-SynCom21-simulated-biomass-in-soil.pdf');
saveas(figG, FigName,'pdf');