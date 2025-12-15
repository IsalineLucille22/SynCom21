function [Sim_Vec, Sim_String_Name, Sorted_ind, Sum_dist] = fun_Similarities_Species(names_species, index_species, Corr_Mat, Euclid_Mat, Hamming_Mat, Cos_sim)
%Function used in ResourceProfile.m. Based on the consumption profils of
%each species, create a consumption similarity index.
Norm_Corr_Mat = (1 - Corr_Mat(index_species,:))/sum(1 - Corr_Mat(index_species,:));
Norm_Euclid_Mat = Euclid_Mat(index_species,:)/sum(Euclid_Mat(index_species,:));
Norm_Hamming_Mat = Hamming_Mat(index_species,:)/sum(Hamming_Mat(index_species,:));
Norm_Cos_sim = (1 - Cos_sim(index_species,:))/sum(1 - Cos_sim(index_species,:));
Sum_dist = Norm_Corr_Mat + Norm_Euclid_Mat + Norm_Hamming_Mat + Norm_Cos_sim;
[Sim_Vec, Sorted_ind] = sort(Sum_dist);
Sim_String_Name = names_species(Sorted_ind);
if Sorted_ind(1) ~= index_species %If the most similar species is not itself
    ind_1 = Sorted_ind(1); ind_2 = Sorted_ind(2);
    Sorted_ind(1) = ind_2; Sorted_ind(2) = ind_1;
end
end