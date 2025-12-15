function lh = logLikelihood_2(Evol_mass, gamma, Time_step, Evol_mass_R, mu_max, rho, T_lag, Ks, opts_1)
tspan = [0 max(Time_step)];
if mu_max>0 && rho>0 && Ks>0 
    sol = ode45(@(t,x) dfun_2(t, x, rho, mu_max, T_lag, Ks), tspan, [Evol_mass(1) Evol_mass_R(1)], opts_1);
    x_est = deval(sol, Time_step);
    x_est = x_est(1,:);
%     temp = -1/(2*gamma^2)*sum((Evol_mass - x_est').^2);
%     lh = sum(temp);  
    lh = -1/2*sum(sum((Evol_mass - x_est').^2./(gamma'.^2),2));
else
lh=-1e30;
end
end