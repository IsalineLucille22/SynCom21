function [Fin_mat, Sigma_mat, h_mat, p_mat, Slope_mat, Intercept_mat] = Val_Trace(Trace_mat, Burn_out_val)
[n, m, l] = size(Trace_mat);
Fin_mat = zeros(n,m);
Sigma_mat = zeros(n,m); h_mat = ones(n,m); p_mat = zeros(n,m);
Slope_mat = ones(n,m); Intercept_mat = zeros(n,m);
maxLags = l - Burn_out_val + 1;
for i = 1:n
    for j = 1:m
        val = reshape(Trace_mat(i,j,Burn_out_val:end), [], 1);
        if std(val) > 0
            % lower_bound = prctile(val, 5); % 5th percentile
            % upper_bound = prctile(val, 95); % 95th percentile
            % val = val(val > lower_bound & val < upper_bound);
            maxLags = length(val) - 1;
            pdf = fitdist(val, 'Normal');
            mu = pdf.mu;
            sigma = pdf.sigma;
            Fin_mat(i,j) = mu;
            Sigma_mat(i,j) = sigma;
            figure(1)
            autocorr(val, 'NumLags', round(l*0.2)); 
            a = polyfit((1:length(val))/1000, val, 1);
            [h_mat(i,j), p_mat(i,j)] = kstest((val - mean(val)) / std(val));
            Slope_mat(i,j) = a(1); Intercept_mat(i,j) = a(2);
            figure(2)
            plot(val)
            hold on
            mean_val = mean(val);
            fplot(@(x) mean_val, [0 maxLags], 'g--')
            figure(3)
            plot(reshape(Trace_mat(i,j,:), 1, []))
            hold on
            mean_val = mean(reshape(Trace_mat(i,j,:), 1, []));
            fplot(@(x) mean_val, [0 l], 'g--')
            % if h_mat(i,j) == 0
            %     c = 10;
            % end
        end
    end
end
end