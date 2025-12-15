function [Mat_Average, diff_obs, diff_rand] = Average_Mat(Struct)
n = length(Struct) - 1;
sum_mat = Struct{1};
[n_row, n_col] = size(sum_mat);
Mat_Average = zeros(n_row, n_col);
for i = 2:n
    sum_mat = sum_mat + Struct{i};
end
for i = 1:n_row
    for j = 1:n_col
        sample = zeros(n,1);
        for k = 1:n
            mat_temp = Struct{k};
            sample(k) = mat_temp(i, j);
        end
        dist_temp = fitdist(sample, 'Normal');
        Mat_Average(i,j) = max(normrnd(dist_temp.mu, dist_temp.sigma),0);%max(random(dist_temp),0);%exprnd(dist_temp.mu);%
    end
end
Mat_Average = 1/n*sum_mat;
end