clear
close all

Name_file = 'data_to_save_liquid_3.mat';
load(strcat('Data/', Name_file));

data_liquid_no_CF = data_to_save_liquid_3;
S = 21;
[data_liquid_no_CF.CrossFeed_Mat_Temp] = deal(zeros(S, S));  % Distributes the zero matrix to all structs


% Threshold_Pred_init = max(normrnd(1.0e-04, 0*0.1e-05, 21, 1),0);
% Pred_Mat_Lyso_init = zeros(S,S);
% Prey_num = [3 7 11 13 15 17 20];
% ind_Lyso = 6;
% Pred_Mat_Lyso_init(ind_Lyso, Prey_num) = 1e-03; %Predation rate randomly chosen
% [data_liquid_rand.Pred_Mat_Lyso] = deal(Pred_Mat_Lyso_init);
% [data_liquid_rand.Threshold_Pred] = deal(Threshold_Pred_init);

save(strcat(cd, '/Data/', 'data_liquid_no_CF'), 'data_liquid_no_CF')