function Delta = cliffsDelta(x, y)
nx = numel(x);
ny = numel(y);

n_sup = 0;
n_inf = 0;
for i = 1:nx
    n_sup = n_sup + sum(y < x(i));
    n_inf = n_inf + sum(y > x(i));
end

Delta = (n_sup - n_inf)/(nx*ny);
end