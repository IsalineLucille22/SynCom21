function dz_vect = fun_CF_Death_Pred(t, z, kappa, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix, Threshold, Threshold_death, Threshold_Pred, Death_Mat, Pred_Mat, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, Hill_CF, Hill_Pred)
z(z<0) = 0;
r = z((end - nb_Res + 1):end);
z = z(1:(end - nb_Res));
z = reshape(z',3, S); 
z = z';%[S_i, P_i, W_i]
k = Hill_CF;%10;
T = (1./(1 + (Threshold(1)./z(:,3)).^k));%(1./(1 + (Threshold./z(:,3)).^k));
T(isnan(T)) = 0;
k = Hill_Pred;%10;
T_death = (1./(1 + (Threshold_death(1)./z(:,1)).^k));%(1./(1 + (Threshold_death(1)./(max(z(:,3) - sum(r),0))).^k));%Threshold on resource concentration?
T_death(isnan(T_death)) = 0;
T_pred = (1./(1 + (Threshold_Pred(1)./z(:,1)).^10));
T_pred(isnan(T_pred)) = 0;
K_S_res = (kappa(:,2) + kappa(:,3))./kappa(:,1).*Resource_Matrix;
K_S_res(K_S_res==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
K_S_Waste = (CrossFeed_Mat + Mat_kappa_3)./kappa(:,1); %Division by columns
K_S_Waste(K_S_Waste==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
K_S_death = (Death_Mat + kappa(:,3))./kappa(:,1);
K_S_death(K_S_death==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
K_S_Pred = (Pred_Mat + kappa(:,3))./kappa(:,1);%(kappa(:,2) + kappa(:,3))./kappa(:,1).*Pred_Mat;%Normal %Comment
K_S_Pred(K_S_Pred==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
dx_vect =  z(:,1).*sum(kappa(:,2).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)),2)...
    + z(:,1).*sum(CrossFeed_Mat.*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)),2)...
    - z(:,1).*sum(Death_Mat.*(T_death.*z(:,3))'./(z(:,3)' + K_S_death), 2)...
    + yield_Pred.*z(:,1).*sum(((t > Lag_time_Pred).*Pred_Mat).*(T_pred.*z(:,1))'./(z(:,1)' + K_S_Pred), 2)...
    - sum(z(:,1).*((t > Lag_time_Pred).*Pred_Mat).*(T_pred.*z(:,1))'./(z(:,1)' + K_S_Pred))';
dy_vect = zeros(S,1); %Turn it into the inhibition byproducts.
dw_vect = z(:,1).*sum(kappa(:,3).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)),2)...
    + z(:,1).*sum(Mat_kappa_3.*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)),2)...
    - sum(z(:,1).*(CrossFeed_Mat + Mat_kappa_3).*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)))'...
    + sum(z(:,1).*Death_Mat.*(T_death.*z(:,3))'./(z(:,3)' + K_S_death))'...
    + (1 - yield_Pred).*z(:,1).*sum(((t > Lag_time_Pred).*Pred_Mat).*(T_pred.*z(:,1))'./(z(:,1)' + K_S_Pred), 2);
dr = - sum(z(:,1).*(kappa(:,2) + kappa(:,3)).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)));
dz_vect = [dx_vect dy_vect dw_vect];
dz_vect = reshape(dz_vect', 1, []);
dz_vect = [dz_vect dr]';
end