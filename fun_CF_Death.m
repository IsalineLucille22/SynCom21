function dz_vect = fun_CF_Death(t, z, kappa, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix, Threshold_CF, Threshold_death, Death_Mat, death_rate, S, Lag_time_Cons, Lag_time_Pred, nb_Res, Hill_CF, Hill_Pred)
%Lag_time_Pred = Lag_time_self_inhibition = 0
z(z<0) = 0;
r = z((end - nb_Res + 1):end);
z = z(1:(end - nb_Res));
z = reshape(z',3, S); 
z = z';%[S_i, P_i, W_i]
k = Hill_CF;
T = (1./(1 + (Threshold_CF./z(:,3)).^k));%(1./(1 + (Threshold_CF(1)./z(:,3)).^k));%
T(isnan(T)) = 0;
k = Hill_Pred;
T_death = (1./(1 + (Threshold_death./z(:,1)).^k));%(1./(1 + (Threshold_death(1)./z(:,1)).^k));
T_death(isnan(T_death)) = 0;
T_long_Term_death = 1/(1 + (0.004./(max(0.0075 - sum(r),0))).^2);
T_long_Term_death(isnan(T_long_Term_death)) = 0;
K_S_res = (kappa(:,2) + kappa(:,3))./kappa(:,1).*Resource_Matrix;
K_S_res(K_S_res == 0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
K_S_Waste = (CrossFeed_Mat + Mat_kappa_3)./kappa(:,1); %Division by columns
K_S_Waste(K_S_Waste == 0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
K_S_death = (Death_Mat + kappa(:,3))./kappa(:,1);
K_S_death(K_S_death == 0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
dx_vect =  z(:,1).*sum(kappa(:,2).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)), 2)...
    + z(:,1).*sum(CrossFeed_Mat.*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)),2)...
    - z(:,1).*sum(((t > Lag_time_Pred).*Death_Mat).*(T_death.*z(:,3))'./(z(:,3)' + K_S_death), 2) ... 
    - z(:,1).*(T_long_Term_death.*death_rate);%z(:,1).*(T_death.*death_rate);%(t > 168).* %Remove 0 multiplication %6.2855e-04 %Modify with time
dy_vect = z(:,1).*sum(CrossFeed_Mat.*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)),2);%zeros(S,1); %Turn it into the inhibition byproducts.
dw_vect = z(:,1).*sum(kappa(:,3).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)),2)...
    + z(:,1).*sum(Mat_kappa_3.*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)),2)...
    - sum(z(:,1).*(CrossFeed_Mat + Mat_kappa_3).*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)))'...
    + sum(z(:,1).*((t > Lag_time_Pred).*Death_Mat).*(T_death.*z(:,3))'./(z(:,3)' + K_S_death))';
dr = - sum(z(:,1).*(kappa(:,2) + kappa(:,3)).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)), 1); %To test
dz_vect = [dx_vect dy_vect dw_vect];
dz_vect = reshape(dz_vect', 1, []);
dz_vect = [dz_vect dr]';
end