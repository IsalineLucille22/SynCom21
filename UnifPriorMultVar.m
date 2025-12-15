function L = UnifPriorMultVar(x, a, b)
x_vect = reshape(x,1,[]);
L = prod(unifrnd(x_vect, a, b));
end