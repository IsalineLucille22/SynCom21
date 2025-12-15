%Calculate specific yield from byproducts per strain and timepoint from the kappa-3 and the CF-matrices

%start by running the simulation script, e.g., Plot_syncom-growth-plots-obs8-model-data_save-w-corr2-start
% change the time step to e.g., 

%Simulations with fitted interspecific interactions
%use here the data_to_save_with_corrections_2 with 10 replicate simulations
%input is S20_S21_abs_abund_cfu_Senka.xlsx, Sheet 4, with 8 replicates 
%starting abundances for Bra, Phe, Mes, Cau, Coh, Tar, Bur, Chi and Muc are corrected

clear
close all

addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data')
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil')

%Name of the data saved
name_Exp = 'Senka_21_CF_Liquid';
Name_file = 'data_to_save_liquid_3.mat';
Name_file_Saved = 'Liquid_System_fitted';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');%%readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8, 'Range','48:69', 'Format','auto');%


Data_Evol = readtable(strcat('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/CoCulture_Soil/Data/Liquid_Data/','abund_cfu_se.xlsx'), 'Sheet', 3, 'Range','26:47', 'Format','auto');%Data in liquid from Clara's experiments
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
Time_step = [0 12 24 48 168]; %Measured time step in hours
S = height(Data_Evol);
Time_step_BP =  0:50:1500;
name = string(table2array(Parameters_set(1:21,1)));
name(end - 1:end) = [{'Ps1'}, {'Ps2'}];

%Load parameters
load(strcat('Data/', Name_file));
data_to_save = data_to_save_liquid_3;
nb_data_set = size(data_to_save,2);

% Fitted parameters
rep = repmat(struct('exp_abs_yield_from_bp', [], 'exp_perc_yield_from_bp', []), 1, nb_data_set);
StackPlotTot = zeros(S, length(Time_step));
StackPlot_Meas_added_errors = zeros(S, length(Time_step)); %We perturbate the observed values accordingling to the simulated variance between data
z = zeros(S, size(data_to_save,2));
z_obs_added_errors = zeros(S, size(data_to_save,2)); %We perturbate the observed values accordingling to the simulated variance between data


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

Data_Evol_temp = table2array(Data_Evol(:, 2:end));
std_y_0 = 0;
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_time_step_BP = length(Time_step_BP);
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
        plot(Time_step, Measured_Abund(j,:,i), '-*', 'Color', colors{j});
        hold on
    end
end
axis([0 180 0 4e-04])
num_fig = num_fig + 1;
rand = unifrnd(0,1);
mean_y_0 = mean(Measured_Abund(:,1,:), 3);
var_mat = var(Measured_Abund, 0, 3);
Measured_Abund_average = mean(Measured_Abund,3);

%Setting for Matlab ODE solver
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.
for zz = 1:nb_data_set 

    ind_sim = zz; 
    kappa_mat = data_to_save(ind_sim).kappa_mat;
    CrossFeed_Mat = data_to_save(ind_sim).CrossFeed_Mat_Temp;
    Resource_Matrix = data_to_save(ind_sim).Resource_Matrix;
    Death_Mat_Temp = data_to_save(ind_sim).Death_Mat_Temp; 
    Threshold_CF = data_to_save(ind_sim).Threshold_CF;
    Threshold_death = data_to_save(ind_sim).Threshold_death; %Change it
    Lag_time_Cons = data_to_save(ind_sim).Lag_time_Cons; 
    Lag_time_Pred = data_to_save(ind_sim).Lag_time_Pred;
    R = data_to_save(ind_sim).R;
    nb_Res = length(R); %Number of resources (group of resources)
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
    death_rate =  data_to_save(ind_sim).death_rate;
    %From Lysobacter
    Pred_Mat_Lyso = data_to_save(ind_sim).Pred_Mat_Lyso;
    Threshold_Pred = data_to_save(ind_sim).Threshold_Pred;  
    
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
    for i = 1:nb_replicates
        %Initial concentrations using a normal distribution
        mat_y_0 = mean_y_0;
        
        mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];
    
        for k = 1:n_exp
            y_0 = sum(mat_y_0(:,1:2),2);
            mat_y_0 = reshape(mat_y_0', 1, []);
            mat_y_0 = [mat_y_0 R];

            sol = ode45(@(t, y) fun_CF_Death_Pred(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix,...
                Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Temp, death_rate,...
                Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
            z_temp = deval(sol, Time_step);
            sum(z_temp(:,1))
            sum(z_temp(:,4))
            sum(z_temp(:,end))
            X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
            P = z_temp(mod(1:S*3,3) == 2,:); %Byproduct growth biomass
            W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
            R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass

            z_temp_BP = deval(sol, Time_step_BP);
            X_BP = z_temp_BP(mod(1:S*3,3) == 1,:); %Species'biomass
            P_BP = z_temp_BP(mod(1:S*3,3) == 2,:); %Byproduct growth biomass
            W_BP = z_temp_BP(mod(1:S*3,3) == 0,:); %Byproducts'biomass
            R_temp_BP = z_temp_BP(S*3 + 1: end,:); %Byproducts'biomass
            Measured_Abund_disp = Measured_Abund_average + normrnd(zeros(S, length(Time_step)), sqrt(var_mat));
    
            z_temp = z_temp(1:(end-nb_Res), end);
            z_temp = reshape(z_temp',3, S);
            z_temp = z_temp';
            z(:,zz) = X(:, end);%Total biomass of all species after 1 week
            z_obs_added_errors(:,zz) = Measured_Abund_disp(:, end);
            StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
            [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
            [Shann_obs, Simp_obs] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));
        end
        StackPlotTot = StackPlotTot + X./sum(X);
        StackPlot_Meas_added_errors = StackPlot_Meas_added_errors + Measured_Abund_disp./sum(Measured_Abund_disp);
   
    end
    
    %Absolute contribution (estimated)
    K_S_Waste = (CrossFeed_Mat + Mat_kappa_3)./kappa_mat(:,1); %Division by columns
    K_S_Waste(K_S_Waste==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
    for i = 1:nb_time_step_BP - 1
        temp = zeros(S);
        for j = 1:S
            T_temp = (1./(1 + (Threshold_CF(1)./W_BP(:,i)).^10));
            %Proportion, for each species, of growth induced by byproducts
            %consumption
            temp(j,:) = (P_BP(j,i + 1) - P_BP(j,i)).*...
                (X_BP(j, i)*CrossFeed_Mat(j,:).*((T_temp.*W_BP(:,i))'./(W_BP(:,i)' + K_S_Waste(j,:))))...
            /sum((X_BP(j, i)*CrossFeed_Mat(j,:).*((T_temp.*W_BP(:,i))'./(W_BP(:,i)' + K_S_Waste(j,:))))); %We can remove X_BP(j, i), because it both appears in denominator and numerato. Here just to relate it to the ODE
        end    
        rep(zz).exp_abs_yield_from_bp{i} = temp;  
    end
    
    %Percentage contribution
    for i = 1:nb_time_step
        temp_2 = zeros(S);  % Initialize the SxS contribution matrix
        for j = 1:S  % Loop over consumers
            T_temp = (1./(1 + (Threshold_CF(1)./W_BP(:,i)).^10));
            BP_contrib = X(j, i)*CrossFeed_Mat(j,:).*((T_temp.*W(:,i))'./(W(:,i)' + K_S_Waste(j,:))); %We can remove X_BP(j, i), because it both appears in denominator and numerato. Here just to relate it to the ODE
            denom = sum(BP_contrib);
            if denom == 0
                denom = 1e-10;  % Small epsilon to avoid NaN
            end
            temp_2(j,:) = BP_contrib/denom;
        end
    
        rep(zz).exp_perc_yield_from_bp{i} = temp_2;
    end
end

%%

close all

[sorted_name, idx] = sort(name); %to sort alphabetically

% take the mean across all simulations for each timepoint

net_bp_yield = cell(1, nb_time_step_BP - 2);
for i = 1:nb_time_step_BP - 2 
    temp = zeros(S, S, nb_data_set);
    
    for j = 1:nb_data_set %simulations
        temp(:,:,j) = rep(j).exp_abs_yield_from_bp{i};    
    end
    
    net_bp_yield{i} = nanmean(temp,3);
end


%make a plot
strains3 = cellfun(@(x) x(1:min(3, end)), sorted_name, 'UniformOutput', false)';

close all

n = 3;%4; % number of rows
m = 3;%7; % number of columns

figH = figure;
asy_ind = zeros(1, n*m);
for i = 2:(n*m + 1)
    % subplot('Position', positions(i,:))
    subplot(n, m, i - 1)
    
    %normalize the display to the maximum value at each time point
    m_u = net_bp_yield{i};
    
    %m_u = mean_usage(:,(i-1)*20+1:(i-1)*20+20)./max(max(mean_usage(:,(i-1)*20+1:(i-1)*20+20)));
    
    % sort to strain alphabetical order
    m_u = m_u(:,idx);
    m_u = m_u(idx,:);
    asy_ind(i - 1) = AsymmetryIndex(m_u);
    
    h = heatmap(strains3, strains3, m_u);
    h.MissingDataColor = [1 1 1]; %white
    h.FontSize = 6;
    h.Title = append('Time',num2str(Time_step_BP(i)));
end
mean_asy_ind = mean(asy_ind);

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'mean-byproduct-strain-yields-per-time.pdf');
saveas(figH,FigName,'pdf');

%plot for a single simulation
%chose a single simulation to plot

%e.g., 

j = 3; %Choice for the parameter sample

close all

figH = figure;
%Total time points value is 28 for abs_yield
n = 3;
m = 3;
for i = 2:(n*m + 1)%nb_time_step
    %subplot('Position',positions(i,:))
    subplot(n, m, i - 1)
    
    m_u = rep(j).exp_abs_yield_from_bp{i};
    
    %sort to strain alphabetical order
    
    m_u = m_u(:,idx);
    m_u = m_u(idx,:);
    
    h = heatmap(strains3, strains3, m_u);
    h.MissingDataColor = [1 1 1]; %white
    h.FontSize = 6;
    h.Title = append('Time',num2str(Time_step_BP(i)));
end

iFolderName = strcat(cd, '/Figures/');
FigName = strcat(iFolderName, name_Exp, 'sim3-byproduct-strain-yields-per-time.pdf');
saveas(figH,FigName,'pdf');