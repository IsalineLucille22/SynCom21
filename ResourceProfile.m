clear;
close all;
clc;

% Similarity between resource preferences based on mono-cultures

%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file = 'Heatmap_resource_preferences';
fig_num = 1;

%Loadind data
Data = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 11,'Range','25:46','Format','auto');%readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 11,'Range','1:21','Format','auto');%
Parameters_set = table2array(Data(:,2:(end - 1))); %Be careful to consider all the indices or not to have some last empty columns
%Parameters_set = Parameters_set/max(max(Parameters_set));
MaxVal = max(max(Parameters_set));%
S = length(Parameters_set(:,1));

figure(fig_num)
Resource_Heatmap = heatmap(Parameters_set);
caxis([0 MaxVal]);
title('Heatmap resource perferences')
fig_num = fig_num + 1;

figure(fig_num)
lowestValue = min(Parameters_set(Parameters_set(:)>0));
highestValue = max(Parameters_set(:));
imagesc(Parameters_set);
cmap = jet(256);
colormap(cmap);
caxis(gca,[lowestValue-2/256, highestValue]);
% Make less than lowest value black:
cmap(1,:)=[.7 .7 .7];
colormap(cmap)
colorbar
yticklabels = table2cell(Data(:,1));
yticks = linspace(1, length(Parameters_set(:,1)), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels(:))
title('Heatmap resource perferences')
fig_num = fig_num + 1;

[Corr_Mat, Euclid_Mat, Hamming_Mat, Cos_sim] = fun_Covariance(Parameters_set./sum(Parameters_set,2)); %Addition of the cosine similarity. Make a dendrogram
figure(fig_num)
imagesc(Corr_Mat);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels(:))
set(gca, 'XTick', yticks, 'XTickLabel', yticklabels(:))
colormap(cmap)
colorbar
title('Correlation matrix')
fig_num = fig_num + 1;
figure(fig_num)
imagesc(Euclid_Mat);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels(:))
set(gca, 'XTick', yticks, 'XTickLabel', yticklabels(:))
colormap(cmap)
colorbar
title('Normalized Euclidean matrix')
fig_num = fig_num + 1;
figure(fig_num)
imagesc(Hamming_Mat);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels(:))
set(gca, 'XTick', yticks, 'XTickLabel', yticklabels(:))
colormap(cmap)
colorbar
title('Hamming matrix on binary matrix')
fig_num = fig_num + 1;
figure(fig_num)
imagesc(Cos_sim);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels(:))
set(gca, 'XTick', yticks, 'XTickLabel', yticklabels(:))
colormap(cmap)
colorbar
title('Cosine similarity')
fig_num = fig_num + 1;

Sim_String_Name_Mat = strings(S);
Sim_Mat = zeros(S);
for i = 1:S
    [Sim_Vec, Sim_String_Name, Sorted_ind, Sum_dist] = fun_Similarities_Species(table2array(Data(:,1)), i, Corr_Mat, Euclid_Mat, Hamming_Mat, Cos_sim);
    Sim_String_Name_Mat(:,i) =  string(Sim_String_Name);
    Sim_Mat(i,:) = Sim_Vec;
end
Sim_String_Name_Mat = Sim_String_Name_Mat';

% figure(fig_num)
% Z = zscore(Parameters_set); %Standardization
% [coeff,score,latent] = pca(Parameters_set(:,1:20));
% Xcentered = score*coeff';
% biplot(coeff(:,1:2),'scores',score(:,1:2),'Color','r','Marker','o');
% fig_num = fig_num + 1;

% evaluation = evalclusters(Parameters_set./sum(Parameters_set,2),"kmeans","silhouette","KList",3:20);
% figure(fig_num)
% plot(evaluation)
% fig_num = fig_num + 1;

Tot_Dist = 0*(1 - Corr_Mat) + 0*(1 - Cos_sim) + Euclid_Mat;
Tot_Dist(1:1+size(Tot_Dist,1):end) = 0;
Tot_Dist = (Tot_Dist + Tot_Dist')/2;
tree = linkage(Parameters_set./sum(Parameters_set, 2),'average', @distfunTest);%'euclidean');%
% leafOrder = optimalleaforder(tree, Tot_Dist);

figure(fig_num)
labels = yticklabels;%flip(Sim_String_Name);
% labels = labels(leafOrder);
dendrogram(tree, 'Labels', labels);
fig_num = fig_num + 1;

figure(fig_num)
Graph_Nut_Profile = graph(Tot_Dist);

% Example threshold value (replace with your desired threshold)
threshold_value = 0.325; % Replace with your threshold value

% Find edges to remove
edges_to_remove = find(Graph_Nut_Profile.Edges.Weight > threshold_value);

% Create a filtered graph by removing edges
filtered_G = rmedge(Graph_Nut_Profile, edges_to_remove);

% Plot the filtered node-edge graph
plot(filtered_G, 'Layout', 'force', 'NodeLabel', labels, 'EdgeLabel', filtered_G.Edges.Weight);

% Customize the appearance of the filtered graph (optional)
title('Filtered Node-Edge Graph');

if save_data == 1
    FolderName = strcat(cd, '/Figures/');
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = num2str(get(FigHandle, 'Number'));
        FigName = strcat('Fig', FigName);
        FigName = strcat(FolderName, FigName, Name_file);
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
    end
    save(strcat(cd, '/Data/', Name_file,'_Corr_Mat.mat'), 'Corr_Mat');
    save(strcat(cd, '/Data/', Name_file,'_Euclid_Mat.mat'), 'Euclid_Mat');
    save(strcat(cd, '/Data/', Name_file,'_Hamming_Mat.mat'), 'Hamming_Mat');
    save(strcat(cd, '/Data/', Name_file,'_Resource_Map.mat'), 'Parameters_set');
end

%% Load interspecific matrices
%Comparison between different fitted sets
%Save or Not
%CONTINUE FROM HERE 
%Load parameters
load(strcat('Data/', Name_file));
data_to_save = data_to_save_Inter_v2;%data_to_save_4;
% Fitted parameters
ind_sim_1 = 2; ind_sim_2 = 2;
%First data set
kappa_mat_1 = data_to_save(ind_sim_1).kappa_mat;
CrossFeed_Mat_1 = data_to_save(ind_sim_1).CrossFeed_Mat_Temp;
Resource_Matrix_1 = data_to_save(ind_sim_1).Resource_Matrix;
Pred_Mat_1 = data_to_save(ind_sim_1).Death_Mat_Temp; 
Threshold_1 = data_to_save(ind_sim_1).Threshold_CF;
Threshold_Pred_1 = data_to_save(ind_sim_1).Threshold_death;
Lag_time_Cons_1 = data_to_save(ind_sim_1).Lag_time_Cons; 
Lag_time_Pred_1 = data_to_save(ind_sim_1).Lag_time_Pred;
R_1 = data_to_save(ind_sim_1).R;


%Second data set
kappa_mat_2 = data_to_save(ind_sim_2).kappa_mat;
CrossFeed_Mat_2 = data_to_save(ind_sim_2).CrossFeed_Mat_Temp;
Resource_Matrix_2 = data_to_save(ind_sim_2).Resource_Matrix;
Pred_Mat_2 = data_to_save(ind_sim_2).Death_Mat_Temp; 
Threshold_2 = data_to_save(ind_sim_2).Threshold_CF;
Threshold_Pred_2 = data_to_save(ind_sim_2).Threshold_death;
Lag_time_Cons_2 = data_to_save(ind_sim_2).Lag_time_Cons; 
Lag_time_Pred_2 = data_to_save(ind_sim_2).Lag_time_Pred;
R_2 = data_to_save(ind_sim_2).R;
%% Graphs for interspecific interactions

close all 

fig_num = 1;

Parameters_set = CrossFeed_Mat_1;

[Corr_Mat, Euclid_Mat, Hamming_Mat, Cos_sim] = fun_Covariance(Parameters_set./sum(Parameters_set,2)); %Addition of the cosine similarity. Make a dendrogram

% evaluation = evalclusters(Parameters_set./sum(Parameters_set,2),"kmeans","silhouette","KList",3:20);
% figure(fig_num)
% plot(evaluation)
% fig_num = fig_num + 1;

tree = linkage(Parameters_set./sum(Parameters_set,2),'average', @distfunTest);%linkage(Parameters_set./sum(Parameters_set,2),'average', @distfunTest);

Tot_Dist = 0*(1 - Corr_Mat) + (1 - Cos_sim) + Euclid_Mat;
Tot_Dist(1:1+size(Tot_Dist,1):end) = 0;
Tot_Dist = (Tot_Dist + Tot_Dist')/2;
%leafOrder = optimalleaforder(tree, Tot_Dist);

figure(fig_num)
d = dendrogram(tree, 'Labels', labels);
fig_num = fig_num + 1;

figure(fig_num)
Graph_Nut_Profile = graph(Tot_Dist);

% Example threshold value (replace with your desired threshold)
threshold_value = 0.5;%0.5; % Replace with your threshold value

% Find edges to remove
edges_to_remove = find(Graph_Nut_Profile.Edges.Weight > threshold_value);
filtered_G = rmedge(Graph_Nut_Profile, edges_to_remove);
plot(filtered_G, 'Layout', 'force', 'NodeLabel', labels);%, 'EdgeLabel', filtered_G.Edges.Weight);
title('Filtered Node-Edge Graph');
fig_num = fig_num + 1;

if save_data == 1
    FolderName = strcat(cd, '/Figures/');
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = num2str(get(FigHandle, 'Number'));
        FigName = strcat('Fig', FigName);
        FigName = strcat(FolderName, FigName, Name_file);
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
    end
end

%% Type of interactions
Parameters_set_Pos = kappa_mat_1(:,2).*CrossFeed_Mat_1;
Parameters_set_neg = Pred_Mat_1;

Precentage_pos = (sum(sum(Parameters_set_Pos>0)) - S)/(S - 1)^2;
Precentage_neg = (sum(sum(Parameters_set_neg>0)) - S)/(S - 1)^2;

figure(fig_num)
boxplot([reshape(Parameters_set_Pos, [], 1) reshape(Parameters_set_neg, [], 1)], 'labels', {'Cross-feeding', 'Predation'})
fig_num = fig_num + 1;
figure(fig_num)
boxplot([reshape(Parameters_set_Pos, [], 1) reshape(Parameters_set_neg, [], 1)], 'labels', {'Cross-feeding', 'Predation'})
ylim([0 2.5])
fig_num = fig_num + 1;
figure(fig_num)
boxplot(Parameters_set_Pos', 'labels', yticklabels) %Per consumer 
ylim([0 10])
fig_num = fig_num + 1;
figure(fig_num)
boxplot(Parameters_set_neg', 'labels', yticklabels) %Per predator
ylim([0 10])
fig_num = fig_num + 1;

if save_data == 1
    FolderName = strcat(cd, '/Figures/');
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = num2str(get(FigHandle, 'Number'));
        FigName = strcat('Fig', FigName);
        FigName = strcat(FolderName, FigName, Name_file);
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
    end
    save(strcat(cd, '/Data/', Name_file,'_Corr_Mat.mat'), 'Corr_Mat');
    save(strcat(cd, '/Data/', Name_file,'_Euclid_Mat.mat'), 'Euclid_Mat');
    save(strcat(cd, '/Data/', Name_file,'_Hamming_Mat.mat'), 'Hamming_Mat');
    save(strcat(cd, '/Data/', Name_file,'_Resource_Map.mat'), 'Parameters_set');
end


%% Similarities between two fitted samples for ineterspecific interactions
close all

nb_species_Comp = S;
Parameters_set = [Pred_Mat_1(1:nb_species_Comp,:);Pred_Mat_2(1:nb_species_Comp,:)];

[Sample_rate, LN_fit, nb_neg_val] = Rates_Distribution(3, {Pred_Mat_1,Pred_Mat_2} );

labels = [yticklabels;yticklabels];

[Corr_Mat, Euclid_Mat, Hamming_Mat, Cos_sim] = fun_Covariance(Parameters_set./sum(Parameters_set,2)); %Addition of the cosine similarity. Make a dendrogram

% evaluation = evalclusters(Parameters_set./sum(Parameters_set,2),"kmeans","silhouette","KList",3:40);
% figure(fig_num)
% plot(evaluation)
% fig_num = fig_num + 1;

tree = linkage(Parameters_set./sum(Parameters_set,2),'average', @distfunTest);%linkage(Parameters_set./sum(Parameters_set,2),'average', @distfunTest);

Tot_Dist = 0*(1 - Corr_Mat) + (1 - Cos_sim) + Euclid_Mat;
Tot_Dist(1:1+size(Tot_Dist,1):end) = 0;
Tot_Dist = (Tot_Dist + Tot_Dist')/2;
leafOrder = optimalleaforder(tree, Tot_Dist);

figure(fig_num)
dendrogram(tree, 'Labels', labels)
fig_num = fig_num + 1;

figure(fig_num)
Graph_Nut_Profile = graph(Tot_Dist);

% Example threshold value (replace with your desired threshold)
threshold_value = 0.5;%0.5; % Replace with your threshold value

% Find edges to remove
edges_to_remove = find(Graph_Nut_Profile.Edges.Weight > threshold_value);
filtered_G = rmedge(Graph_Nut_Profile, edges_to_remove);
h_filtered = plot(filtered_G, 'Layout', 'force', 'NodeLabel', labels);%, 'EdgeLabel', filtered_G.Edges.Weight);
title('Filtered Node-Edge Graph');
fig_num = fig_num + 1;

Sim_String_Name_Mat_CF = strings(2*nb_species_Comp);
nb_sim = 0;
nb_sim_1 = 0;
nb_no_sim = 0;
for i = 1:(2*nb_species_Comp)
    [Sim_Vec_CF, Sim_String_Name_CF, Sorted_ind_CF] = fun_Similarities_Species(labels, i, Corr_Mat, Euclid_Mat, Hamming_Mat, Cos_sim);
    Sim_String_Name_Mat_CF(:,i) =  string(Sim_String_Name_CF);
    if Sim_String_Name_Mat_CF(1,i) == Sim_String_Name_Mat_CF(2,i)
        nb_sim = nb_sim + 1;
    elseif sum(contains(Sim_String_Name_Mat_CF(2:9,i), Sim_String_Name_Mat_CF(1,i))) > 0
        nb_sim_1 = nb_sim_1 + 1;
    elseif sum(contains(Sim_String_Name_Mat_CF(21:end,i), Sim_String_Name_Mat_CF(1,i))) == 0
        nb_no_sim = nb_no_sim + 1;
    end
end
Sim_String_Name_Mat_CF = Sim_String_Name_Mat_CF';
nb_sim 
nb_sim_1
nb_no_sim