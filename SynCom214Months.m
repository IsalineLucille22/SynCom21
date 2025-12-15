%Script made by Isaline Guex
clear;
close all;

%Save or Not
save_data = 0; %1 if save, 0 otherwise

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','94:114', 'Format','auto');
Data_Evol = readtable(strcat('Data/SynCom21_prolonged_incubation/Up_to_4mon_incubation_SynCom21/','Niche_Flow_cfu_tax_SynCom21.xlsx'), 'Sheet', 1, 'Range','1:841', 'Format','auto');
Rep_name = ["Syn3m_1", "Syn3m_2", "Syn3m_3", "Syn3m_4"];
len_Rep_name = length(Rep_name);
S = height(unique(Data_Evol.SynComGenus));
 nb_fig = 1;

slope = zeros(S - 1, len_Rep_name);
for k = 1:len_Rep_name
    Data_Eval = table2array(Data_Evol(Data_Evol.Replicate == Rep_name(k),4:8)); %Change 4:8
    Time_step = unique(Data_Eval(:,2)*24); %Time in hours
    name_species = table2array(unique(Data_Evol(:,2)));
    nb_Time_step = length(Time_step);
    %Loadind data
    Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46', 'Format','auto');
    Sorted_name = string(table2array(Parameters_set(1:S-1,1)));%string(table2array(table(2:22,3)));
    
    %Initialization of the colors
    colors = distinguishable_colors(60);
    
    Data = zeros(S, nb_Time_step);
    
    figure(nb_fig)
    for i = 1:S
        Data(i,:) = Data_Eval(((i - 1)*nb_Time_step + 1): i*nb_Time_step, 5);
        plot(Time_step, Data(i,:), '-*', 'Color', colors(i,:))
        hold on
    end
    
    nb_fig = nb_fig + 1;
    
    Data(Data == 0) = 1e-10;
    Data = Data(:,2:nb_Time_step);
    log_data = log(Data);
    
    ind_Lyso = name_species == "Lysobacter";
    name_species(ind_Lyso) = [];
    log_data(ind_Lyso,:) = [];
    for i = 1:(S - 1)
        p = polyfit(Time_step(2:nb_Time_step), log_data(i,:)',1);
        order_temp = find(string(name_species(i)) == Sorted_name);
        slope(order_temp, k) = p(1);
    end
end


slope(slope > 0) = -1e-05;