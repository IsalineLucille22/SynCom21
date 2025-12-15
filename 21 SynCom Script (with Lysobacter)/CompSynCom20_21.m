%Script to compare abundances of the SynCom20 member with or without
%Lysobacter
clear;
close all;

% %Save or Not
% save_data = 0; %1 if save, 0 otherwise

%%Loading data
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','142:163', 'Format','auto');
Data_Evol_SynCom21 = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto');
Data_Evol_SynCom20 = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
Time_step = [0 1 3 7 10 21]*24;%[0 12 48 96 7*24 21*24];%Measured time step in hours
Time_step_BP = 1:50:500; %Definition of a time step for the byproduct visualization
tspan = [0, max(Time_step)]; %Time interval in hours
S = height(Data_Evol_SynCom21);
Data_Evol_SynCom21_wo_Lysobacter = Data_Evol_SynCom21(1:S ~= 6, :);
S = S - 1;
num_fig = 1;

Data_Evol_temp_1 = table2array(Data_Evol_SynCom21_wo_Lysobacter(:, 2:end));% - table2array(Data_Evol_SynCom20(:, 2:end));
nb_obs = length(Data_Evol_temp_1(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
Measured_Abund_1 = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund_1(:,:,i) = Data_Evol_temp_1(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
figure(num_fig)
StackPlot_Meas_1 = Measured_Abund_1./sum(Measured_Abund_1);
StackPlot_Meas_mean_1 = mean(StackPlot_Meas_1, 3);
bar(StackPlot_Meas_mean_1(:, end))
num_fig = num_fig + 1;

Data_Evol_temp_2 = table2array(Data_Evol_SynCom20(:, 2:end));
nb_obs = length(Data_Evol_temp_2(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
Measured_Abund_2 = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund_2(:,:,i) = Data_Evol_temp_2(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
figure(num_fig)
StackPlot_Meas_2 = Measured_Abund_2./sum(Measured_Abund_2);
StackPlot_Meas_mean_2 = mean(StackPlot_Meas_2, 3);
bar(StackPlot_Meas_mean_2(:, end))
num_fig = num_fig + 1;

figure(num_fig)
diff = StackPlot_Meas_mean_2 - StackPlot_Meas_mean_1;%>0 higher SynCom20, preys > 0
bar(diff(:, end))

num_fig = num_fig + 1;
figure(num_fig)
diff_2 = Measured_Abund_2 - Measured_Abund_1; diff_2 = mean(diff_2, 3);%>0 higher SynCom20, preys > 0
bar(diff_2(:, end))
num_fig = num_fig + 1;

figure(num_fig)
bar(diff)
xlabel('Row index (n)')
ylabel('Value')
legend('T0','T1','T2', 'T3', 'T4', 'T5')
title('Grouped Bar Plot of Matrix Rows')