function [New_Mat] = Increased_Mat(Mat, ind, n_row, n_col, ind_square)
New_Mat = zeros(n_row, n_col);
if ind_square == 1
    New_Mat = zeros(n_row, n_col);
    New_Mat([1:(ind - 1) (ind + 1):end],1:(ind - 1)) = Mat(1:end,1:(ind - 1));
    New_Mat(1:(ind - 1), (ind + 1):end) = Mat(1:(ind - 1),ind:end);
    New_Mat((ind + 1):end,(ind + 1):end) = Mat(ind:end, ind:end);
else
    New_Mat(1:(ind - 1),:) = Mat(1:(ind - 1),:);
    New_Mat((ind + 1):end,:) = Mat(ind:end,:);
end
end