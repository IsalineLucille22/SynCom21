clear;
close all;

%Save or Not
save_data = 1; %1 if save, 0 otherwise
Name_file = 'Death_Model_Loop_2';

%Loadind data
Data = readtable(strcat('Data/','ResourcesDepletionClara.xlsx'), 'Sheet', 1,'Range','1:485', 'Format','auto');
name_Resources = Data(:,1);
nb_res = size(name_Resources);
nb_res = nb_res(1)/4;
Res_Conc = table2array(Data(:,3:4));

colors = distinguishable_colors(nb_res + 1);
tot_res = zeros(4,1);

figure(1)
for i = 1:nb_res
    temp = Res_Conc(4*(i-1) + 1: 4*i, :);
    plot(temp(:, 1), temp(:, 2), '-o', 'Color', colors(i,:));
    hold on
    tot_res = tot_res + temp(:, 2);
end

figure(2)
plot(temp(:, 1), tot_res, '-o', 'Color', colors(i + 1,:));