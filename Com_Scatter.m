function num_fig = Com_Scatter(title_name, z_fin_obs, z_fin_sim, num_fig, Select_S, colors, name, lim_x, lim_y)
figure(num_fig);
S = length(z_fin_obs); S_temp = length(Select_S); name_temp = name(Select_S);
z_fin_obs_temp = z_fin_obs(Select_S); z_fin_sim_temp = z_fin_sim(Select_S); 
colors_temp = colors(Select_S, :); colors_temp_comp = colors(S+1-Select_S,:);
for j = 1:S_temp
    scatter(z_fin_obs_temp(j), z_fin_sim_temp(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors_temp(j,:), 'MarkerFaceColor', colors_temp_comp(j,:));     
    hold on
end
axis([0 lim_x 0 lim_y]);
reflin = refline(1,0);
axis square
reflin.Color = 'r';
xlabel('Experiment'); 
ylabel('Simulation');
legend(name_temp);
title(title_name)
num_fig = num_fig + 1;
end