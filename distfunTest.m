function D2 = distfunTest(ZI,ZJ)
m2 = length(ZJ(:,1));
D2 = zeros(m2,1);
for j = 1:m2
%         D2(j) = sqrt(sum((ZI/sum(ZI) - ZJ(j,:)/sum(ZJ(j,:))).^2));
        D2(j) = sqrt(sum((ZI - ZJ(j,:)).^2)) + 0*(1 - ZI*ZJ(j,:)'/(norm(ZI)*norm(ZJ(j,:))));%(1 - Cos_sim) + Euclid_Mat
        %Cos_sim(i,j) = Mat_Res_Pref(i,:)*Mat_Res_Pref(j,:)'/(norm(Mat_Res_Pref(i,:))*norm(Mat_Res_Pref(j,:))); % Cosine similarity
end
end