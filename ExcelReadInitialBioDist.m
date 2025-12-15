% Read Excel files and generate figures showing the impact of a 100-fold increase 
% in one species (while keeping total biomass constant). Compare results across 
% 10 parameter sets obtained using the Metropolis-Hastings (MH) algorithm.
clear
close all


%Save or Not
save_data = 0; %1 if save, 0 otherwise

%Loadind data
Data_set = readtable(strcat('Data/','InitialBiomassDisturbances.xlsx'), 'Sheet', 2,'Range','1:21', 'Format','auto');
S = height(Data_set);
Data = table2array(Data_set(:, 2:12));
name = table2array(Data_set(:,1));
values_heterogeneity = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; %0 = standard system, every species could interact with the other. 1 = no CF. For 0 and 1 there is no stochasticity

figure(1)
bar(name, Data);
legend(string(values_heterogeneity), 'Orientation', 'vertical', 'Location', 'southeast')

figure(2)
bar(values_heterogeneity', Data', 'stacked');
xlabel('Cross-feeding heterogeneity')
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
