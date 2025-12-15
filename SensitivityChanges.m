function New_Mat = SensitivityChanges(Old_Mat, cov_mat)
[n, m] = size(Old_Mat);
New_Mat = Old_Mat + normrnd(zeros(n, m), cov_mat);%mvnrnd(meanVector, covarianceMatrix, numSamples);
% mean_vec = zeros(m,1);
% to_add = mvnrnd(mean_vec, cov_mat, n);
% New_Mat = Old_Mat + to_add;
New_Mat = max(New_Mat, 0);
end