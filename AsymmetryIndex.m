function ass_ind = AsymmetryIndex(mat)
 diff = mat - mat';
 ass_ind = norm(diff, 'fro')/norm(mat, 'fro');
end