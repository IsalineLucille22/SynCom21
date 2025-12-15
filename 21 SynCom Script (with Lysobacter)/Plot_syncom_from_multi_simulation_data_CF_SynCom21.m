%Simulations with fitted interspecific interactions
%use here the data_to_save_no_CF with 4 replicate simulations
%input is S20_S21_abs_abund_cfu_Senka.xlsx, Sheet 4, with 8 replicates 
%starting abundances for Bra, Phe, Mes, Cau, Coh, Tar, Bur, Chi and Muc are corrected


%CHANGE INITIAL ABUNDANCES BRUNA'S DATA


clear
close all

% addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
% addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')
addpath('/Users/pret_helpdesk/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
addpath('/Users/pret_helpdesk/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')


name_Exp = 'Bruna_21_CF_Reduce_Res_Init_Senka';
Name_file = 'struct_tot_SynCom21_Soil_only_Lyso_Parfor';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Newv4';%'struct_tot_SynCom21_Soil_only_Lyso_Parfor';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Newv3';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_New';%'data_Reduced_dim';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Newv2';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_New';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_New';%'data_to_save_SynCom21_Soil_only_Lyso_Parfor_Random';%'data_to_save_Inter_v2';%'data_no_CF.mat';
Name_file_Resources_Death = 'data_to_save_SynCom21_Reduced'; %Data fitted on moncultures

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');

%actual absolute abundance data, has quadriplicate values for each time point

% Data_Evol = readtable(strcat('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/Liquid_Data/','abund_cfu_IG.xlsx'), 'Sheet', 2, 'Range','46:67', 'Format','auto');%1:22 without correction for 0. Data in soil extract from Phil and Clara's experiments
% Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 7, 'Range','1:22', 'Format','auto'); %Senka's data
Data_Evol = readtable(strcat('Data/','SSC21_genera_relative-abundances.xlsx'), 'Sheet', 4, 'Range','1:22', 'Format','auto'); %Bruna's data
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
S = height(Data_Evol);
Time_step = [0 12 22 38 70 168 504];%Bruna %[0 1 3 7 10 21]*24;%Senka %[0 12 48 96 7*24 21*24];%Clara %Measured time step in hours
nb_res = 12;

%Load parameters
load(strcat('Data/', Name_file));
data_to_save = struct_tot_SynCom21_Soil_only_Lyso_Parfor;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_Newv4;%struct_tot_SynCom21_Soil_only_Lyso_Parfor;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_Newv3;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_New;%data_Reduced_dim; %data_to_save_SynCom21_Soil_only_Lyso_Parfor_New;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_New;%data_to_save_SynCom21_Soil_only_Lyso_Parfor_7;%data_to_save_Inter_v2;%data_no_CF;
load(strcat('Data/', Name_file_Resources_Death));
data_to_save_Res_Death = data_to_save_Mono;

% Fitted parameters
StackPlotTot = zeros(S, length(Time_step));
StackPlot_Meas_added_errors = zeros(S, length(Time_step)); %We perturbate the observed values accordingling to the simulated variance between data
z = zeros(S, size(data_to_save,2));
z_obs_added_errors = zeros(S, size(data_to_save,2)); %We perturbate the observed values accordingling to the simulated variance between data
exp_biomass = zeros(S, length(Time_step), size(data_to_save,2));
exp_resources = zeros(nb_res, length(Time_step), size(data_to_save,2));
exp_byproducts = zeros(S, length(Time_step), size(data_to_save,2));


%Initialization of the colors
colors_ind = [4, 1, 10, 15, 12, 3, 8, 6, 20, 18, 19, 2, 5, 21, 13, 9, 14, 7, 16, 17, 11];
colors_init = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#B6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#91665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC', '#4393C3', '#B35806'};%distinguishable_colors(S);
colors = {};
for i = 1:S
    colors{i} = colors_init{colors_ind(i)};
end

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
%display the biomass variations among the 4 observed replicates
num_fig = 1;
figure(num_fig)
for i = 1:nb_rep
    for j = 1:S
        plot(Time_step, Measured_Abund(j,:,i), '-*', 'Color', colors{j});%, Time_step, R_temp, 'o');
        hold on
    end
end
axis([0 600 0 4e-04])
num_fig = num_fig + 1;
rand = unifrnd(0,1);
mean_y_0 = mean(Measured_Abund(:,1,:), 3);
var_mat = var(Measured_Abund, 0, 3);
Measured_Abund_average = mean(Measured_Abund,3);

%Setting for Matlab ODE solver
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-12);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
nb_data_set = size(data_to_save,2);
[Shann_sim, Simp_sim, Shann_obs, Simp_obs] = deal(zeros(nb_data_set, nb_time_step));
for zz = 1:nb_data_set

    ind_sim = zz; 
    kappa_mat = data_to_save(ind_sim).kappa_mat;
    CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
    Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; %data_to_save_Res_Death.Death_Mat_Temp; %data_to_save(ind_sim).Death_Mat_Temp; 
    Threshold_CF = data_to_save(ind_sim).Threshold_CF;
    Threshold_death = data_to_save(ind_sim).Threshold_death; 
    Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; %data_to_save_Res_Death.Lag_time_Cons; %data_to_save(ind_sim).Lag_time_Cons; 
    % Lag_time_Cons_Pseudo = data_to_save_Res_Death.Lag_time_Cons;
    % Lag_time_Cons_Pseudo(20:21,:) = Lag_time_Cons(20:21,:);
    Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
    death_rate =  data_to_save(ind_sim).death_rate;
    name = string(table2array(Parameters_set(1:21,1)));
    Pred_Mat_Lyso = data_to_save(ind_sim).Pred_Mat_Lyso;
    Threshold_Pred = data_to_save(ind_sim).Threshold_Pred;
    R = data_to_save(ind_sim).R/1.5; %Reduce 1.25 + init Senka. 1.5 just reduce
    nb_Res = length(R); %Number of resources (group of resources)
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix; %data_to_save_Res_Death.Resource_Matrix; %data_to_save(ind_sim).Resource_Matrix; %Addition of a line for Lysobacter
    % Resource_Matrix_Pseudo = data_to_save_Res_Death.Resource_Matrix; %data_to_save(ind_sim).Resource_Matrix; %Addition of a line for Lysobacter
    % Resource_Matrix(20:21,:) = Resource_Matrix_Pseudo(20:21,:);
    
    
    % Measured_Abund = table2array(Data_Evol(1:20, 2:7));
    % Number of surviving species after 8 weeks
    Threshold_Surviving = 1e-10;
    nb_Surv_Obs = sum(Measured_Abund_average > Threshold_Surviving);
    
    %Initialization of the model parameters fixed for all replicates 
    t_0 = 0; %Time 0
    n_exp = 1; %No transfer in Philip experiment 
    tspan = [0, max(Time_step)*24]; %Time interval in hours
    yield_Pred = 0.2;% 20% of yield for predation
    
    nb_replicates = 1;
    num_fig = 2;
    rand_rep = randi(nb_rep);
    mean_y_0 = Measured_Abund(:, 1, rand_rep);
    % mean_y_0(16) = 0; %For Bruna's data, no microbacterium
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
            X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
            P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
            W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
            R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass

            figure(num_fig)
            for j = 1:S
                plot(Time_step, X(j,:), '-o', 'Color', colors{j});%, Time_step, R_temp, 'o');
                hold on
            end
            legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
            axis([0 600 0 4e-04]) %axis([0 600 0 2e-03])
            num_fig = num_fig + 1;
            figure(num_fig)
            Measured_Abund_disp = Measured_Abund_average + normrnd(zeros(S, length(Time_step)), sqrt(var_mat));
            for j = 1:S
                plot(Time_step, Measured_Abund_disp(j,:), '-*', 'Color', colors{j});%, Time_step, R_temp, 'o');
                hold on
            end
            legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
            axis([0 600 0 4e-04])
    
            num_fig = num_fig + 1;
            z_temp = z_temp(1:(end-nb_Res), end);
            z_temp = reshape(z_temp', 3, S);
            z_temp = z_temp';
            z(:,zz) = X(:, end);%X(:, 4);%Total biomass of all species after 1 week or at the end of the experiement
            z_obs_added_errors(:,zz) = Measured_Abund_disp(:, end);
            StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
            %Diversity indices
            [Shann_sim(zz, :), Simp_sim(zz, :)] = Shannon_Simpson_Indices(S, StackPlot);
            [Shann_obs(zz, :), Simp_obs(zz, :)] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));
        end
        StackPlotTot = StackPlotTot + X./sum(X);
        StackPlot_Meas_added_errors = StackPlot_Meas_added_errors + Measured_Abund_disp./sum(Measured_Abund_disp);
    
        exp_biomass(:,:, zz) = X;
        exp_resources(:,:, zz) = R_temp; 
        exp_byproducts(:,:, zz) = W;
    end

end
mean_Shann_sim = mean(Shann_sim); mean_Shann_obs = mean(Shann_obs);
mean_Simp_sim = mean(Simp_sim); mean_Simp_obs = mean(Simp_obs);


var_data = var(exp_biomass, 0, 3);

StackPlotTot = StackPlotTot/zz;
StackPlot_Meas_added_errors = StackPlot_Meas_added_errors/zz;

StackPlotTot = StackPlotTot./nb_replicates;
z_fin_sim = mean(z, 2); %Absolute stationary abundances
z_fin_obs_added_errors = mean(z_obs_added_errors, 2);%z_fin_obs = Measured_Abund(1:S, end); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
disp(sum(z_fin_sim)/sum(z_fin_obs_added_errors))

h_vect = zeros(1,S);
p_vect = zeros(1,S);
for i = 1:S
    [h_vect(i),p_vect(i)] = ttest2(StackPlot_Meas_added_errors(i,:)', StackPlotTot(i,:)');
end
name_diff = name(logical(h_vect));

figure(num_fig);
bar(StackPlotTot', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;

figure(num_fig); 
bar(mean(StackPlot_Meas_added_errors, 3)', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
num_fig = num_fig + 1;

[diff_vect, I, fact] = Sort_diff(StackPlotTot(:,end), StackPlot_Meas_added_errors(:,end));
name_order = name(I);
h = figure(num_fig);
stem(diff_vect);
xtickangle(90)
set(gca,'xtick',1:21,'xticklabel',name_order)
ylabel('Sim - Obs')
num_fig = num_fig + 1;

figure(num_fig);
for j = 1:S
    scatter(StackPlot_Meas_added_errors(j, n_exp+1), StackPlotTot(j,n_exp+1), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
    hold on
end
axis([0 0.7 0 0.7]);
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
    scatter(z_fin_obs_added_errors(j), z_fin_sim(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
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

iFolderName = strcat(cd, '/Figures/');
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = num2str(get(FigHandle, 'Number'));
    FigName = strcat('Fig', FigName);
    FigName = strcat(iFolderName, FigName, name_Exp, 'SynCom_21_');
    set(0, 'CurrentFigure', FigHandle);
    saveas(FigHandle, FigName, 'pdf');
end
%% Computation of the rates distribution for each species

% [~, nb_iter] = size(data_to_save_SynCom21_Soil_only_Lyso_Parfor_7);
% Mat_Struct = {data_to_save_SynCom21_Soil_only_Lyso_Parfor_7.CrossFeed_Mat_Temp};
% [Sample_rate, LN_fit, nb_neg_val] = Rates_Distribution(12, Mat_Struct, num_fig, name);

%% plot figure with the observed or simulated biomass per species at the end

%cd('/Users/jvanderm/Library/CloudStorage/OneDrive-UniversitédeLausanne/SynCom model paper/Data')

close all
figK = figure;
subplot('Position',[0.3 0.2 0.62 0.2]);
bar(mean(Measured_Abund(:, end, :), 3))

hold on

for i = 1:nb_rep %number of replicates
    scatter(1:S, Measured_Abund(:, end, i),60,'.','r')
end

set(gca, 'YScale', 'log')
ylim([1e-7 3e-4])
yticks([1e-6 1e-5 1e-4]);

xticks(1:length(name)); 
xticklabels(name);
grid on
ylabel('Biomass g/mL')
title('observed');

subplot('Position',[0.3 0.5 0.62 0.2]);

bar(mean(z,2))
hold on
for i=1:nb_rep %number of replicates
 scatter(1:S, z(:,i), 60, '.', 'k')
end

set(gca, 'YScale', 'log')
ylim([1e-7 3e-4])
grid on
ylabel('Biomass g/mL')
yticks([1e-6 1e-5 1e-4]);
title('simulated');

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'expected-end-biomass-SynCom20-in-soil-data-CF.pdf');
% saveas(figK,FigName,'pdf');
print(figK,FigName, '-dpdf', '-painters');


% sort of double scatter plot for simulated and experimental data
figL = figure;
subplot('Position',[0.1 0.1 0.55 0.55]);

for i=1:nb_rep %number of replicates
    scatter(median(z,2),  Measured_Abund(:, end, i), 60,'.','k')
    hold on
end
axis([0 2e-4 0 2e-4])
axis square
hold on
for i = 1:nb_data_set %strains
    y =  median(Measured_Abund(:, end, :),3);
    scatter(z(:,i), y, 60,'.','b')
    hold on
end

set(gca, 'YScale', 'log', 'XScale','log')
ylim([1e-7 3e-4])
xlim([1e-7 3e-4])
reflin = refline(1,0);
axis square
reflin.Color = 'r';
ylabel('Experiment'); 
xlabel('Simulation');
grid on
xticks([1e-6 1e-5 1e-4]);


iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'expected-end-biomass-SynCom20-in-soil-scatter.pdf');
% saveas(figL,FigName,'pdf');
print(figL,FigName, '-dpdf', '-painters');

%individual boxcharts

close all

figL = figure;
subplot('Position',[0.1 0.1 0.55 0.55]);
boxchart(z')
set(gca, 'YScale', 'log')
ylim([1e-7 3e-4])
%saveas(figL,'boxchart-z.pdf','pdf');


close all

figL = figure;
subplot('Position',[0.1 0.1 0.55 0.55]);
boxchart(reshape(Measured_Abund(:, end, :), S, nb_rep)')
set(gca, 'YScale', 'log')
% ylim([1e-7 3e-4])
%saveas(figL, strcat(iFolderName, name_Exp, 'boxchart-observed.pdf'),'pdf');
print(figL, strcat(iFolderName, name_Exp, 'boxchart-observed.pdf'), '-dpdf', '-painters');

close all

figK = figure;

tiledlayout(3,7)

for i = 1:S  
    nexttile(i)
    
    for k = 1:nb_data_set
        plot(Time_step, exp_biomass(i, :, k),'m-');
        hold on
    end
    hold on
    
    for kk = 1:nb_rep %number of replicate observations
        plot(Time_step, Measured_Abund(i, :, kk),'r--');
    end
    
    set(gca,'YScale','log');
    ylim([0 5e-4]);
    xlim([0 500]);
    yticklabels({}); %turn off the ytick labels, but they go from 1e-8 to 1e-4
    grid on
    bact_name = name(i);
    title(bact_name);
end

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'ind-growth-plots-model-data_CF-start.pdf');
%saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');

%% display mean and stdev of the cross-feeding matrices over time

tmp = zeros(S, S, nb_data_set);

for i = 1:nb_data_set
    tmp (:,:,i) = data_to_save(i).CrossFeed_Mat_Temp;
end

mean_CF = mean(tmp,3);
var_CF = std(tmp,[],3);

close all

figK = figure;

subplot('Position',[0.1 0.1 0.3 0.3]);

imagesc(mean_CF');
colorbar

title('crossfeeding');

subplot('Position',[0.5 0.1 0.3 0.3]);

imagesc(var_CF');
colorbar

title('crossfeeding std');

%saveas(figK,'heatmap-sim-cross-feeding-data_to_save_w_corr.pdf','pdf');
%%Useless if considering SynCom21 because the variations only occur for
%%Lysobacter

% multi stack bar plot of relative abundances in the simulations
% compared to the mean of 8 observations

[strains, idx]=sortrows(name); %to sort alphabetically
mean_CF=mean_CF(idx,:); %sort the byproducts rows accordingly
mean_CF=mean_CF(:,idx); %sort the columns accordingly

var_CF=var_CF(idx,:);
var_CF=var_CF(:,idx);

close all

figK=figure;

subplot('Position',[0.15 0.2 0.3 0.4]);
h = heatmap(strains, strains, mean_CF, Colormap=hot(8));

h.Title = 'mean cross-feeding';

subplot('Position',[0.6 0.2 0.3 0.4]);
h = heatmap(var_CF, Colormap=pink);

h.Title = 'cross-feeding standard error';

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'heatmap-names-sim-cross-feeding-data_to_save_w_corr2.pdf');
% saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');

close all

figK = figure;

%tiledlayout(3,7)
tiledlayout(6,7)

for i = 1:nb_data_set
    nexttile(i)
    b = bar((exp_biomass(:,:,i)./sum(exp_biomass(:,:,i)))','stacked','FaceColor','flat');%b = bar((exp_biomass{i}./sum(exp_biomass{i}))','stacked','FaceColor','flat');
    for j=1:S
        b(j).FaceColor = colors{j};
    end
    title(strcat('simulation',num2str(i)));
end

StackPlot_Meas_m = mean(StackPlot_Meas, 3);

nexttile
b = bar(StackPlot_Meas_m','stacked','FaceColor','flat');
for j=1:S
    b(j).FaceColor = colors{j};
end
title('mean obs');

legend(name, 'Orientation', 'vertical', 'Location', 'none','Position',[0.75 0.1 0.2 0.4])

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'rel-abund-model-data-CF.pdf');
% saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');

%%simple plot of the maximum biomass attained in the simulations versus observations

sim_fin_biomass = zeros(nb_data_set, 1);
for i = 1:nb_data_set
    % sim_max_biomass(i,:) = max(sum(exp_biomass(:,:,i)));
    sim_fin_biomass(i,:) = max(sum(exp_biomass(:,end,i)));
end

%repeat calculation of the measured abundances
obs_fin_biomass = zeros(nb_rep, 1);

for i = 1:nb_rep
    % Obs_Abund = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
    % obs_max_biomass(i,:) = max(sum(Obs_Abund));
    % obs_fin_biomass(i,:) = max(sum(Measured_Abund(:,:,i)));
    obs_fin_biomass(i,:) = max(sum(Measured_Abund(:,end,i)));
end

close all
figH = figure;

subplot('Position',[0.1 0.1 0.2 0.4])

y = [mean(sim_fin_biomass); mean(obs_fin_biomass)];

X = categorical({'Sim','Obs'});
X = reordercats(X,{'Sim','Obs'});
b = bar(X,y,'FaceColor','flat');
set(gca,'FontSize',7);
b.CData(1,:) = [1 0 1];
b.CData(2,:) = [0 1 1];

hold on
for j = 1:nb_data_set
    scatter(X(1), sim_fin_biomass(j),50,'k','.')
end
hold on
for j = 1:nb_rep
 scatter(X(2), obs_fin_biomass(j),50,'k','.')
end
ylabel('biomass (gC/mL)')
ylim([0 2e-3])

title('max biomass','FontSize',9);

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'max-data-CF-sim-obs-biomass.pdf');
% saveas(figH,FigName,'pdf');
print(figH,FigName, '-dpdf', '-painters');

[Sortedname, idx] = sort(name);
z = z(idx, :);
save(strcat('Data/',name_Exp, 'z'), 'z')

%% plot resource usage

close all

for i = 1:nb_data_set
p = plot(Time_step, exp_resources(:,:,i)');
	for j = 1:nb_Res
	    p(j).Color = colors{j};
	end
hold on
end


%% plot byproduct formation by strain for all replicates

close all

figK = figure;

% tiledlayout(4,5)

tiledlayout(3, 7)
% tiledlayout(6, 7)

for i=1:S

    nexttile(i)

	for k = 1:nb_data_set
	    % plot(Time_step, exp_byproducts{k}(i,:),'m-');
        plot(Time_step, exp_byproducts(i,:, k),'m-');
	    hold on
	end
	hold on

    %set(gca,'YScale','log');
    ylim([0 0.004]);
    yticklabels({}); %turn off the ytick labels, but they go from 0 to 0.004
    grid on
    bact_name=name(i);
    title(bact_name);
end

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'ind-byproducts-data_save_w_corr2-start.pdf');
% saveas(figK,FigName,'pdf');
print(figK,FigName, '-dpdf', '-painters');

% as heatmap with the mean of the byproducts at each time point

tmp = exp_byproducts;

% for i = 1:10 %nr of replicates
% 
% 	tmp(:,:,i)= exp_byproducts{i};
% 
% end

mean_byproducts=mean(tmp, 3);

[strains, idx]=sortrows(name); %to sort alphabetically
mean_byproducts=mean_byproducts(idx,:); %sort the byproducts matrix accordingly

close all

figK=figure;

subplot('Position',[0.2 0.2 0.2 0.55]);
h = heatmap(Time_step/24, strains, mean_byproducts);

h.Title = 'Byproduct formation';
h.XLabel = 'Time (days)';


iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'byproducts-heatmap.pdf');
% saveas(figK,FigName,'pdf');
print(figK, FigName, '-dpdf', '-painters');

%% sum byproducts in the system plus the modeled resource utilization

figK = figure;

subplot('Position',[0.1 0.1 0.3 0.3]);

	for k = 1:nb_data_set
	    % byproduct_sum_t = sum(exp_byproducts_t{k},1);
	    % resource_sum_t = sum(exp_resources_t{k},1);
	    byproduct_sum = sum(exp_byproducts(:, :, k),1);
	    resource_sum = sum(exp_resources(:, :, k),1);
    
	    plot(Time_step/24,byproduct_sum,'k-')
    
	    hold on
	    plot(Time_step/24,resource_sum,'m-')
	end

xlabel('Time (days)');
ylabel('Concentration (mg C/ml)');
grid on


iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'sum-byproducts-resources.pdf');
% saveas(figK,FigName,'pdf');
print(figK,FigName, '-dpdf', '-painters');

%% Growth rates plots 

% cd('/Users/jvanderm/Library/CloudStorage/OneDrive-UniversitédeLausanne/SynCom model paper/Data')


% individual growth plots

%plot(Time_step,X)

% obs_X = [2, 6, 10, 14, 18, 22];%[2,10,18,26,34,42,50,58];

% Time_step = [0 1 3 7 10 15 21 22]*24; %for the experimental data


%plot(Time_step,obs_X)

close all

figK = figure;

% tiledlayout(3, 7)
tiledlayout(6, 7)

for i = 1:S

    nexttile(i)

	for k = 1:size(exp_biomass,2)
	    % plot(Time_step, exp_biomass{k}(i,:),'m-');
        plot(Time_step, exp_biomass(i,:,k),'m-');
	    hold on
	end
	hold on

    for kk = 1:nb_rep %number of replicate observations
        plot(Time_step, Measured_Abund(i, :, kk),'r--');
    end

	% for kk = 1:8 %number of replicate observations
	%     plot(Time_step2, table2array(Data_Evol(i, obs_X-1+kk)),'b-');
	% end

    set(gca,'YScale','log');
    ylim([1e-8 5e-4]);
    yticklabels({}); %turn off the ytick labels, but they go from 1e-8 to 1e-4
    grid on
    bact_name=name(i);
    title(bact_name);
end

%%calculate mu growth rates

%%Post processing analysis with moving interval on LN-transformed area data to calculate mu-Monod

mu_sim = zeros(nb_data_set, S);

for i = 1:nb_data_set
    A = exp_biomass(:,:,i)';%exp_biomass{i}';
    
    %%extract mu and rsq with a gliding scale - here 4 points  
    t = Time_step; %time scale in hours for the length of the table 
    int = 4; %number of points for the floating interval
    
    for k = 1:S %size(A,2)
	    for z = 1:length(t) - int
	        [mu_max(z), slope(z), intercept(z)] = regression(t(z:z+int), log(A(z:z+int, k))'); 
	    end
        max_m = max(mu_max);%m(r>0.8)); %take the max value of from all intervals and where r2>0.95
        mu_sim(i,k) = max_m;
        % clearvars best_m   %Doesn't needed. Will be overwrite anyway  
    end
end %end pos loop


mu_exp = zeros(nb_rep, S);
for i = 1:nb_rep %what

    % A = table2array(Data_Evol(:,(obs_X)-1+i))';
    A = Measured_Abund(:, :, i)';
    
    %%extract mu and rsq with a gliding scale - here 4 points
    t = Time_step; %Time_step2; %time scale in hours for the length of the table
    
    int = 3; %number of points for the floating interval
    figure(2)
    for k = 1:S %size(A,2)
	    for z = 1:length(t)-int
            temp_data = A(z:z+int,k);
            temp_data(temp_data <= 0) = 1e-07;
	        [mu_max, slope, intercept] = regression(t(z:z+int), log(temp_data)');%regression(t(z:z+int), log(A(z:z+int,k))'); 
            plot(t(z:z+int), log(temp_data), '-o')
            hold on
	    end
        max_m = max(mu_max);%max(m(r>0.5)); %take the max value of from all intervals and where r2>0.95
        mu_exp(i,k) = max_m;
        %clearvars best_m 
    end
end %end pos loop


%% Growth rates computation 

%addpath('/Users/jvanderm/Documents/Tania Research/microLife_manuscript/Miguel_Trabajo_Zenodo_data/Surface_model/Model/version2/Figures generation/')
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/Shared spatial/SurfaceModels/Figures generation')


%For observations
for i = 1:nb_rep 

    A = Measured_Abund(:, :, i)';
    [time_init, mu_log] = deal(zeros(S, 1));
    for k = 1:S
    
        Evol_mass = A(:,k);
        lag_time = 0;
    
        x0 = [1, max(Evol_mass), Evol_mass(1), lag_time];
        opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-10);
        fun_Logistic = @(x, Time_step) x(2)./(1 + (x(2) - x(3))/x(3).*exp(-x(1)*Time_step));
    
        [pars, ~] = fminsearch(@(x) -logLikelihood(Evol_mass, 1e02, Time_step, 0, x(1), x(4), x(2), opts_1),...
                x0,...%Initial guess
                optimset('MaxFunEvals', 20000));
    
        sol = ode45(@(t,x) dfun(t, x, Evol_mass(1), 0, pars(1), pars(4), pars(2)), [0 max(Time_step)], Evol_mass(1), opts_1);
        x_est = deval(sol, Time_step);
    
        time_init(k) = pars(4);
        mu_log(k) = pars(1);
    
    end     
       pos.time_init{i} = time_init; %fitted lag time (h)
       pos.mu_log{i} = mu_log; %fitted logaritmic growth rate - to compare to mu_monod; multiply by log(2)
end %end fitting loop


%%with mu-Monod?

%For simulations
[mu_max_Monod_LT, LT_Monod] = deal(zeros(nb_data_set, S));
for k = 1:nb_data_set %number of replicate simulations

    data = exp_biomass(:, :, k); 
    for i = 1:S
    
        data_Evol = data(:,1:end);
        Evol_mass = data_Evol(i, :);
        Evol_mass_R = sum(R_temp); %mean_R_0*ones(1, length(Time_step));
        gamma = std(data_Evol); %Standard deviation of the observations. 
        opt_1 = odeset('RelTol',1e-9,'AbsTol',1e-10); %Setting for ODEs solver
                
        %Implicit Monod estimation with lag time
        r_1 = 1 ;
        Ks = max(Evol_mass)/2; %Half velocity constant converted in equivalent OD value 0.0581 for OD.
        x0 = [0.05, 0.1, 0.5, Ks]; %[mu_max, rho, t_lag, Ks]
        pars_2 = fminsearch(@(x) -logLikelihood_2(Evol_mass, gamma, Time_step, Evol_mass_R, x(1), x(2), x(3),  x(4), opt_1),...
            x0,...%Initial guess
            optimset('MaxFunEvals', 20000));
        sol = ode45(@(t,x) dfun_2(t, x, pars_2(2), pars_2(1), pars_2(3),  pars_2(4)), [0 max(Time_step)], [Evol_mass(1) Evol_mass_R(1)], opt_1);
        x_est_2 = deval(sol, Time_step);
        x_est_2 = x_est_2(1,:);
        
        mu_max_Monod_LT(k,i) = pars_2(1); %to have per hour
        LT_Monod(k,i) = pars_2(3); %to have hour    
    end
end

%%remove any mu-prediction above 1

mu_max_Monod_LT(mu_max_Monod_LT > 1)=nan;

[mu_max_Monod_exp, LT_Monod_exp] = deal(zeros(nb_rep, S));
for k = 1:nb_rep %number of experimental repetitions

    % data=table2array(Data_Evol(:,(obs_X)-1+k))';
    data = Measured_Abund(:, :, k)';
    
    for i = 1:S %the strains
    
        data_Evol = data(:, 1:end); %Needed?
        Evol_mass = data_Evol(:,i);
        Evol_mass_R = sum(R_temp);%mean_R_0*ones(1, length(Time_step2)); %We d
        gamma = std(data_Evol,0,2); %Standard deviation of the observations. 
        opt_1 = odeset('RelTol',1e-9,'AbsTol',1e-10); %Setting for ODEs solver
                
        %Implicit Monod estimation with lag time
    
        Ks = max(Evol_mass)/2;% %Half velocity constant converted 
    
        x0 = [0.05, 0.1, 1, Ks]; %[mu_max, rho, t_lag, Ks]
        pars_2 = fminsearch(@(x) -logLikelihood_2(Evol_mass, gamma, Time_step, Evol_mass_R, x(1), x(2), x(3),  x(4), opt_1),...
            x0,...%Initial guess
            optimset('MaxFunEvals', 20000));
        sol = ode45(@(t,x) dfun_2(t, x, pars_2(2), pars_2(1), pars_2(3),  pars_2(4)), [0 max(Time_step)], [Evol_mass(1) Evol_mass_R(1)], opt_1);
        x_est_2 = deval(sol, Time_step);
        x_est_2 = x_est_2(1,:);
        
        mu_max_Monod_exp(k, i) = pars_2(1); %to have per hour
        LT_Monod_exp(k, i) = pars_2(3); %to have hour  
    end
end

%%remove any mu-prediction above 1

mu_max_Monod_exp(mu_max_Monod_exp > 1) = nan;

% %% Questions about this section
% %Are those growth rates the previous ones fitted?
% %plot growth rate predictions in SynCom
% 
% %data to use is the file Growth_rates_SynCOM_soil_indiv_mix_sims.xlsx
% 
% % cd('/Users/jvanderm/Library/CloudStorage/OneDrive-UniversitédeLausanne/SynCom model paper/Data/Data IG')
% cd('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/SynCom model paper/Data/Data IG')
% % clear
% % close all
% 
% %Change accordingly to the SynCom21
% strains = readtable('Growth_rates_SynCOM_soil_indiv_mix_sims.xlsx','Range','A2:A22');
% [strains, idx] = sortrows(strains); %sort alphabetically
% 
% %read the 10 replicate mu-Monod from the simulations
% 
% mu_sim=readtable('Growth_rates_SynCOM_soil_indiv_mix_sims.xlsx','Range','C2:L22');
% mu_sim=mu_sim(idx,:);
% 
% mu_sim=table2array(mu_sim);
% mu_sim(mu_sim>0.35)=nan;
% 
% %read the 8 replicate mu-Monod from the simulations
% 
% mu_exp_obs = readtable('Growth_rates_SynCOM_soil_indiv_mix_sims.xlsx','Range','N2:U22');
% mu_exp_obs = mu_exp_obs(idx,:); %reorganize according to alphabetical strain name
% 
% mu_exp_obs = table2array(mu_exp_obs);
% mu_exp_obs(mu_exp_obs > 0.35) = nan;
% 
% %read mu of individual growth
% 
% mu_ind=readtable('Growth_rates_SynCOM_soil_indiv_mix_sims.xlsx','Range','AB2:AD22');
% mu_ind=mu_ind(idx,:); %reorganize according to alphabetical strain name
% mu_ind=table2array(mu_ind);
% 
% % make a plot
% 
% xvalues=categorical(table2array(strains));
% X=reordercats(xvalues, table2array(strains));
% 
% close all
% figK=figure;
% 
% subplot('Position',[0.3 0.2 0.62 0.4]);
% 
% bar(X,nanmean(mu_sim,2))
% 
% hold on
% 
% for i=1:10 %number of replicates
% scatter([1:20],mu_sim(:,i),60,'.','b')
% end
% 
% 
% bar(X,nanmean(mu_exp,2))
% hold on
% for i=1:8 %number of replicates
% scatter([1:20],mu_exp(:,i),60,'.','g')
% end
% 
% scatter(X,nanmean(mu_ind,2),'+');
% hold on
% for i=1:3 %number of replicates
% scatter([1:20],mu_ind(:,i),60,'.','r')
% end
% 
% 
% %ylim([1e-7 3e-4])
% grid on
% ylabel('max growth rate h-1')
% title('max growth rates');
% legend({'blue=sim','green=obs','cross=individual'},'Location','eastoutside');
% 
% saveas(figK,'obs-sim-mu_max-SynCom20-in-soil-data-w-corr2.pdf','pdf');
% 
% 
% % ttest differences
% 
% mu_test=mattest(mu_sim,mu_exp);
% 
% FDR = mafdr(mu_test, 'BHFDR', true); % Estimate positive false discovery rate for multiple hypothesis testing
% 
% mu_test=[strains table(mu_test) table(FDR)];
% 
% writetable(mu_test, 'mu-sim-obs-mattest.csv');