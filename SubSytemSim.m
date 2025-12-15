%Simulations subsytems with fitted interspecific interactions
clear
close all


%Save or Not
save_data = 0; %1 if save, 0 otherwise

%Loadind data
%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','94:114', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');%readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');%
Data_Evol_test = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
S = height(Data_Evol);


%Load parameters
Time_step = [0 1 3 7 10 21]*24;%0:1:168;%Time step in hours [0 1 3 7 10 21 31 61 92 122];%
Name_file = 'data_to_save';
load(strcat('Data/', Name_file));
data_to_save = data_to_save_4;
ind_sim = 3;
kappa_mat = data_to_save(ind_sim).kappa_mat;
CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
Resource_Matrix = data_to_save(ind_sim).Resource_Matrix;
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
Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);

Threshold_Surviving = 1e-10;
Ratio_fin = zeros(1,20);

nb_iter = 1;
Av_diff_fin_tot_temp_abs = zeros(S, nb_iter); Av_diff_fin_tot_temp = zeros(S, nb_iter);
mean_biomass = zeros(S, nb_time_step,  nb_iter); std_biomass = zeros(S, nb_time_step);
Cluster_props = 0;
for n = 1:nb_iter
    %Initialization of the model parameters fixed for all replicates 
    t_0 = 0; %Time 0
    nb_Res = length(R); %Number of resources (group of resources)
    n_exp = 1; %No transfer in Philip experiment 
    tspan = [0, max(Time_step)]; %[0, 22*24]; %Time interval in hours
    S_sim = 15; %Species present into the subsystem
    S = length(S_sim); %Number of species
    Stat_Biomass_Alone = table2array(Parameters_set(1:20,6));
    Stat_Biomass_Alone = Stat_Biomass_Alone(S_sim);
    Av_diff_fin_tot_temp_abs = zeros(S, nb_iter); Av_diff_fin_tot_temp = zeros(S, nb_iter);
    mean_biomass = zeros(S, nb_time_step,  nb_iter); std_biomass = zeros(S, nb_time_step);
    
    upper_tri = rand(S);
    upper_tri = upper_tri > Cluster_props;
    upper_tri = triu(upper_tri);
    upper_tri = upper_tri + upper_tri' - diag(diag(upper_tri));
    % upper_tri = ones(S, S);
    
    kappa_mat_temp = kappa_mat(S_sim,:);
    CrossFeed_Mat_temp = CrossFeed_Mat(S_sim, S_sim).*upper_tri;
    Mat_kappa_3_temp = Mat_kappa_3(S_sim, S_sim);
    Resource_Matrix_temp = Resource_Matrix(S_sim,:);
    Pred_Mat_temp = Pred_Mat(S_sim,S_sim);
    Threshold_temp = Threshold(S_sim);
    Threshold_Pred_temp = Threshold_Pred(S_sim);
    Lag_time_Cons_temp = Lag_time_Cons(S_sim,:);
    Lag_time_Pred_temp = Lag_time_Pred(S_sim,S_sim);
    name_temp = name(S_sim);
    yield_Pred_temp = 0;% 20% of yield for predation
    mean_y_0_temp = mean_y_0(S_sim);
    temp_bio = mean_y_0(S_sim);
    death_rate_temp = death_rate(S_sim);
    % mean_y_0_temp(1) = 0*mean_y_0_temp(1);
    % temp_bio(16) = 100*temp_bio(16);
    % mean_y_0_temp = temp_bio*sum(mean_y_0_temp)/sum(temp_bio);
    
    %Setting for Matlab ODE solver
    nb_tot_Species = S*3 + nb_Res;
    opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
    
    %Initialization of the colors
    
    colors = distinguishable_colors(60);
    
    nb_replicates = 1;
    
    num_fig = 1;
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0_temp;
    
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];
    
    for k = 1:n_exp
        y_0 = sum(mat_y_0(:,1:2),2);
        mat_y_0 = reshape(mat_y_0', 1, []);
        mat_y_0 = [mat_y_0 R];
    
        sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat_temp, CrossFeed_Mat_temp, Mat_kappa_3_temp, Resource_Matrix_temp, Threshold_temp, Threshold_Pred_temp, Pred_Mat_temp, death_rate_temp, S, Lag_time_Cons_temp, Lag_time_Pred_temp, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); 
        z_temp = deval(sol, Time_step);
        % sum(z_temp(:,1))
        % sum(z_temp(:,4))
        % sum(z_temp(:,end))
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
        figure(num_fig)
        for j = 1:S
            plot(Time_step, X(j,:), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
            hold on
        end
        legend(name_temp, 'Orientation', 'vertical', 'Location', 'southeast')
        num_fig = num_fig + 1;
        
        figure(num_fig)
        plot(Time_step, sum(R_temp), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
        num_fig = num_fig + 1;

        figure(num_fig)
        for j = 1:S
            plot(Time_step, W(j,:), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
            hold on
        end
        num_fig = num_fig + 1;
    
    
        z_temp = z_temp(1:(end-nb_Res), end);
        z_temp = reshape(z_temp',3, S);
        z_temp = z_temp';
        z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
        StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
        [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
    end
    StackPlotTot = X./sum(X);
    
    StackPlotTot = StackPlotTot./nb_replicates;
    z_fin_sim = z(1:S,end); %Absolute stationary abundances
    mean_biomass(:, :, n) = X;
    nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
    Ratio_Biomass_Soil_Obs = Stat_Biomass_Alone(2)./z_fin_sim(2);%Test it
    z_fin_obs = mean(Measured_Abund(1:S, end, :), 3); %Absolute stationary abundances
    diff_fin = z_fin_obs - z_fin_sim;
    diff_fin_abs = abs(z_fin_obs - z_fin_sim);
    Av_diff_fin_tot_temp_abs(:, n) = diff_fin_abs;
    Av_diff_fin_tot_temp(:, n) = diff_fin;
    
    
    % figure(num_fig);
    % bar(StackPlotTot', 'stacked');
    % axis([0 11.5 0 1])
    % legend(name_temp, 'Orientation', 'vertical', 'Location', 'southeast')
    % title('Stacked all replicates')
    % num_fig = num_fig + 1;
    % 
    % 
    % figure(num_fig);
    % plot(1:length(Time_step), Simp_sim, 'r--o')
    % num_fig = num_fig + 1;
    % 
    % figure(num_fig);
    % plot(1:length(Time_step), nb_Surv_Sim, 'r--o')
    % num_fig = num_fig + 1;
    % Ratio_fin(n) = Ratio_Biomass_Soil_Obs;
    
    % sum(X)
end
std_biomass = std(mean_biomass, [], 3);
mean_biomass = mean(mean_biomass, 3);
figure(num_fig)
for j = 1:S
    errorbar(Time_step, mean_biomass(j,:), std_biomass(j,:), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
    hold on
end
axis([0 600 0 3e-04])
legend(name_temp, 'Orientation', 'vertical', 'Location', 'southeast')
num_fig = num_fig + 1;

std_diff_fin_tot = std(Av_diff_fin_tot_temp_abs, [], 2);
Av_diff_fin_tot =  sum(Av_diff_fin_tot_temp_abs, 2)/nb_iter;
Av_diff_fin_tot_temp = reshape(Av_diff_fin_tot_temp, [], 1);
figure(num_fig)
histfit(Av_diff_fin_tot_temp)
axis([0.25*min(Av_diff_fin_tot_temp) 1.2*max(Av_diff_fin_tot_temp) 0 1200])
fitdist = fitdist(Av_diff_fin_tot_temp, 'normal');
mu = fitdist.mu; std = fitdist.sigma;
text(min(Av_diff_fin_tot_temp)*0.2, 400, ['(mu, std): (' num2str(mu) ', ' num2str(std) ')'], 'Color','red','FontSize',12)

iFolderName = strcat(cd, '/Figures/');
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = num2str(get(FigHandle, 'Number'));
    FigName = strcat('Fig', FigName);
    FigName = strcat(Name_file, FigName, 'Impact_cluster', num2str(Cluster_props*100));
    set(0, 'CurrentFigure', FigHandle);
    saveas(FigHandle, FigName, 'pdf');
end