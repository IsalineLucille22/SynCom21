function [diff_vect, I, fact] = Sort_diff(Sim, obs)
diff = Sim - obs;
fact = obs./Sim;
[diff_vect, I] = sort(diff);
end