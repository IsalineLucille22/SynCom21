%Script made by Adline Vouillamoz and Isaline Guex
clear;
close all;

%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file_Lyso = 'data_to_save_SynCom21_Soil_only_Lyso_Parfor_7';%'data_to_save_SynCom21_Soil_with_pred_Updated_Parfor_1';%'data_to_save_SynCom21_Soil_with_pred_Updated_Parfor';%

%%Loading data
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','142:163', 'Format','auto');
Data_Evol = readtable(strcat('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/Liquid_Data/','abund_cfu_IG.xlsx'), 'Sheet', 2, 'Range','46:67', 'Format','auto');%1:22 without correction for 0. Data in soil extract from Phil and Clara's experiments
%Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto');
Time_step = [0 12 48 96 7*24 21*24];%[0 1 3 7 10 21]*24;%Measured time step in hours
Time_step_BP = 1:50:1500; %Definition of a time step for the byproduct visualization
tspan = [0, max(Time_step_BP)]; %Time interval in hours
S = height(Data_Evol);

%%%%Classification species growth according to Melanie 
%%%Highly abundant
High_Abund = [6 11 15 16 17 20 21];
%%%Lowly abundant
Low_Abund = [3 7 8 9 10 19];
%%%Very lowly abundant
Very_Low_Abund = [1 2 4 5 12 13 14 18];

    
%Load parameters
load(strcat('Data/', Name_file_Lyso));
data_to_save = data_to_save_SynCom21_Soil_only_Lyso_Parfor_7;%data_to_save_SynCom21_Soil_with_pred_Updated_Parfor_1;%data_to_save_SynCom21_Soil_with_pred_Updated_Parfor;%
% Fitted parameters
ind_sim = 1; %1 to 6 data sets possible, 20 not so good, 1 good
kappa_mat = data_to_save(ind_sim).kappa_mat;
CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; 
Threshold_CF = data_to_save(ind_sim).Threshold_CF;
Threshold_death = data_to_save(ind_sim).Threshold_death; 
Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
death_rate =  data_to_save(ind_sim).death_rate;
name = string(table2array(Parameters_set(1:21,1)));
Pred_Mat_Lyso = data_to_save(ind_sim).Pred_Mat_Lyso;
Threshold_Pred = data_to_save(ind_sim).Threshold_Pred;
R = data_to_save(ind_sim).R; 
nb_Res = length(R); %Number of resources (group of resources)
Resource_Matrix = data_to_save(ind_sim).Resource_Matrix; %Addition of a line for Lysobacter

% %%%%Simulation with the mean of multiple datasets
% values = {data_to_save_SynCom21_Soil_only_Lyso_Parfor_5.CrossFeed_Mat_Temp}; % Collect the field values into a cell array
% tot_mat = values{1, 1};
% for i = 2:length(values)
%     tot_mat = tot_mat + values{1, i};
% end
% CrossFeed_Mat =  tot_mat/length(values);

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
StackPlot_Meas = mean(StackPlot_Meas, 3);
rand = unifrnd(0,1);
mean_y_0 = mean(Measured_Abund(:,1,:), 3);
% mean_y_0 = Measured_Abund(:, 1, 2) +  (1 - rand)*Measured_Abund(:, 1, 1);
% mean_y_0 = mean_y_0 + normrnd(0, 0*0.05*mean(mean_y_0));
Measured_Abund = mean(Measured_Abund, 3);


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

        sol = ode45(@(t, y) fun_CF_Death_Lyso(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix,...
            Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate,...
            Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
        z_temp = deval(sol, Time_step);
        % sum(z_temp(:,1))
        % sum(z_temp(:,4))
        % sum(z_temp(:,end))
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Byproduct growth biomass
        W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
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
        legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 600 0 4e-04])
        num_fig = num_fig + 1;
        figure(num_fig)
        for j = 1:S
            plot(Time_step, Measured_Abund(j,:), '--*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
            hold on
        end
        legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 600 0 4e-04])
        num_fig = num_fig + 1;
        figure(num_fig)
        for j = 1:S
            plot(Time_step_BP, X_BP(j,:), '--*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
            hold on
        end
        legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 2000 0 4e-04])

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


% %To use these graph, we have to change the ODE dy_vect to produce the
% %byproducts. But this ODE has to be equal to 0 when fitting the parameters.
% K_S_Waste = (CrossFeed_Mat + Mat_kappa_3)./kappa_mat(:,1); %Division by columns
% K_S_Waste(K_S_Waste==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
% % %%LOOK AT THESE GRAPHS
% for i = 2:(length(Time_step_BP) - 1)
%     temp = zeros(S);
%     for j = 1:S
%         temp(j,:) = (P_BP(j,i + 1) - P_BP(j,i)).*(X_BP(j, i)*CrossFeed_Mat(j,:).*(W_BP(:,i)'./(W_BP(:,i)' + K_S_Waste(j,:))))/sum((X_BP(j, i)*CrossFeed_Mat(j,:).*(W_BP(:,i)'./(W_BP(:,i)' + K_S_Waste(j,:)))));
%     end
%     figure(num_fig)
%     num_fig = num_fig + 1;
%     heatmap(temp)
% end
% 
% 
% %Percentage contribution
% for i = 1:length(Time_step) 
%     temp_2 = zeros(S);  % Initialize the SxS contribution matrix
%     for j = 1:S  % Loop over consumers
%         BP_contrib = X(j, i)*CrossFeed_Mat(j,:).*(W(:,i)'./(W(:,i)' + K_S_Waste(j,:)));
%         denom = sum(BP_contrib);
%         if denom == 0
%             denom = 1e-10;  % Small epsilon to avoid NaN
%         end
%         temp_2(j,:) = BP_contrib/denom;
%     end
%     figure(num_fig)
%     num_fig = num_fig + 1;
%     heatmap(temp_2)
% end



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

num_fig = Com_Scatter('High Abund Relative', StackPlot_Meas(:, end), StackPlotTot(:, end), num_fig, High_Abund, colors, name, 1.01*max([StackPlot_Meas(High_Abund, end); StackPlotTot(High_Abund, end)]), 1.01*max([StackPlot_Meas(High_Abund, end); StackPlotTot(High_Abund, end)]));
num_fig = Com_Scatter('Low Abund Relative', StackPlot_Meas(:, end), StackPlotTot(:, end), num_fig, Low_Abund, colors, name, 1.01*max([StackPlot_Meas(Low_Abund, end); StackPlotTot(Low_Abund, end)]), 1.01*max([StackPlot_Meas(Low_Abund, end); StackPlotTot(Low_Abund, end)]));
num_fig = Com_Scatter('Very Low Abund Relative',StackPlot_Meas(:, end), StackPlotTot(:, end), num_fig, Very_Low_Abund, colors, name, 1.01*max([StackPlot_Meas(Very_Low_Abund, end); StackPlotTot(Very_Low_Abund, end)]), 1.01*max([StackPlot_Meas(Very_Low_Abund, end); StackPlotTot(Very_Low_Abund, end)]));

num_fig = Com_Scatter('High Abund Absolute', z_fin_obs(:, end), z_fin_sim(:, end), num_fig, High_Abund, colors, name, 1.01*max([z_fin_obs(High_Abund, end); z_fin_sim(High_Abund, end)]), 1.01*max([z_fin_obs(High_Abund, end); z_fin_sim(High_Abund, end)]));
num_fig = Com_Scatter('Low Abund Absolute',z_fin_obs(:, end), z_fin_sim(:, end), num_fig, Low_Abund, colors, name, 1.01*max([z_fin_obs(Low_Abund, end); z_fin_sim(Low_Abund, end)]), 1.01*max([z_fin_obs(Low_Abund, end); z_fin_sim(Low_Abund, end)]));
num_fig = Com_Scatter('Very Low Abund Absolute',z_fin_obs(:, end), z_fin_sim(:, end), num_fig, Very_Low_Abund, colors, name, 1.01*max([z_fin_obs(Very_Low_Abund, end); z_fin_sim(Very_Low_Abund, end)]), 1.01*max([z_fin_obs(Very_Low_Abund, end); z_fin_sim(Very_Low_Abund, end)]));

figure(num_fig);
plot(1:length(Time_step), Simp_obs, 'b--o')
hold on
plot(1:length(Time_step), Simp_sim, 'r--o')
num_fig = num_fig + 1;

iFolderName = strcat(cd, '/Figures/');
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = num2str(get(FigHandle, 'Number'));
    FigName = strcat('Fig', FigName);
    FigName = strcat(iFolderName, FigName, 'SynCom_21');
    set(0, 'CurrentFigure', FigHandle);
    saveas(FigHandle, FigName, 'pdf');
end

X = X'; %For R sript comparison