function [Corr_Mat, Euclid_Mat, Hamming_Mat, Cos_sim] = fun_Covariance(Mat_Res_Pref)
%Function used in ResourceProfile.m. Computation of different metrics to
%assess the similarity between species consumption profil.
[S, R] = size(Mat_Res_Pref);
Binary_Mat = zeros(S,R);
Binary_Mat(Mat_Res_Pref > 0) = 1;
% Cov_Mat = zeros(S);
% for i = 1:S
%     for j = i:S
% %         Cov_Mat(i,j) = 1/(R-1)*sum((Mat_Res_Pref(i,:) - mean(Mat_Res_Pref(i,:))).*(Mat_Res_Pref(j,:) - mean(Mat_Res_Pref(j,:))))/sqrt(var(Mat_Res_Pref(i,:))*var(Mat_Res_Pref(j,:)));
%         Cov_Mat(i,j) = 1/(R-1)*(Mat_Res_Pref(i,:) - mean(Mat_Res_Pref(i,:)))*(Mat_Res_Pref(j,:) - mean(Mat_Res_Pref(j,:)))'/sqrt(var(Mat_Res_Pref(i,:))*var(Mat_Res_Pref(j,:)));
%     end
% end
% Cov_Mat = Cov_Mat + Cov_Mat' - eye(S);
%Correlation matrix
Corr_Mat = corr(Mat_Res_Pref');

%Euclidean distance
Euclid_Mat = zeros(S);
Cos_sim = zeros(S);
for i = 1:S
    for j = i:S
%         Euclid_Mat(i,j) = sqrt(sum((Mat_Res_Pref(i,:)/sum(Mat_Res_Pref(i,:)) - Mat_Res_Pref(j,:)/sum(Mat_Res_Pref(j,:))).^2));
        Euclid_Mat(i,j) = sqrt(sum((Mat_Res_Pref(i,:) - Mat_Res_Pref(j,:)).^2));
%         dist = norm(Mat_Res_Pref(i,:) - Mat_Res_Pref(j,:));
        Cos_sim(i,j) = Mat_Res_Pref(i,:)*Mat_Res_Pref(j,:)'/(norm(Mat_Res_Pref(i,:))*norm(Mat_Res_Pref(j,:))); % Cosine similarity
    end
end
Euclid_Mat = Euclid_Mat + Euclid_Mat';
Cos_sim = Cos_sim + Cos_sim'- eye(S);
Hamming_Mat =  pdist(Binary_Mat, 'hamming');
Hamming_Mat = squareform(Hamming_Mat);
end