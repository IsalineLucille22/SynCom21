function dz = dfun_2(t, x, rho, mu_max, T_lag, Ks)
% Monod function %Replace by concentration?
x(x<=0) = 0;
alpha = 1;% - 1/(1 + (T_alpha/t)^k_alpha);
R_conc = alpha*x(2); %The concentration correspond to x(2), resource biomass which is in g/mL
% KS = mu_max/(rho*Ks);
dz_1 = (t >= T_lag)*x(1)*mu_max*R_conc/(R_conc + Ks); %Include lag time %(kappa_2 + kappa_3)./kappa_1
dR = -(t >= T_lag)*x(1)*(mu_max/rho)*R_conc/(R_conc + Ks);
dz = [dz_1; dR];
end