function drho_vect = Emergent_Encoding_fun(t, rho_vect, x_vect, death_rate, lambda)
rho_vect(rho_vect<0) = 0;
drho_vect = -death_rate.*rho_vect + lambda.*(1 - rho_vect)*x_vect;
end