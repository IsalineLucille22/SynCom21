function out_struct = Weighted_Struct(struct_input)
[~, nb_obs] = size(struct_input);
weights = randi(11, 1, nb_obs) - 1;
weights = weights./sum(weights); %Weights 
fields_names = fieldnames(struct_input);
% nb_fields = length(fields_names);
out_struct = struct();
for i = 1:23
    temp_name = fields_names{i};
    temp_average = weights(1)*struct_input(1).(fields_names{i});
    for j = 2:nb_obs
        temp_average = temp_average + weights(j)*struct_input(j).(fields_names{i});
    end
    out_struct.(temp_name) = temp_average;
end
end