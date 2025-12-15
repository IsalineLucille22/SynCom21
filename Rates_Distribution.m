function [Sample_rate, LN_fit, nb_neg_val] = Rates_Distribution(ind_species, Mat_Struct, num_fig, name)
%Generate the histrogram and log-normal fit for the rates of the
%interaction matrices of Mat_Struct.
N = length(Mat_Struct);
mat_0 = Mat_Struct{1};
m = length(mat_0(1,:));
Sample_rate = zeros(1, N*m);
name_species = name(ind_species);
for i = 1:N
    mat_temp = Mat_Struct{i};
    Sample_rate((i - 1)*m + 1: i*m) = mat_temp(ind_species,:);
end
non_neg_val = Sample_rate(Sample_rate>0);
nb_neg_val = N*m - length(non_neg_val);
figure(num_fig)
LN_fit = lognfit(non_neg_val);
histfit(non_neg_val, [], 'lognormal')
axis([0 20 0 10]);
histogram(Sample_rate, 50)

% Generate values for the fitted log-normal curve
x_values = linspace(min(non_neg_val), max(non_neg_val), 100); % Define x-axis values
pdf_values = lognpdf(x_values, LN_fit(1), LN_fit(2)); % Compute the probability density function

% Scale the PDF to match the histogram
histogram_values = histcounts(Sample_rate(Sample_rate > 0), 'Normalization', 'pdf');
scale_factor = max(histogram_values) / max(pdf_values); % Scale the PDF to match the histogram height
pdf_values_scaled = pdf_values * scale_factor;

% Plot the histogram
figure(num_fig);
histogram(Sample_rate, 50, 'Normalization', 'pdf'); % Normalize histogram to plot as PDF
hold on;

% Plot the fitted log-normal curve
plot(x_values, pdf_values_scaled, 'r-', 'LineWidth', 2); % Red curve for log-normal fit
hold off;

% Customize the axes and labels
xlabel('Sample Rate');
ylabel('Density');
title(strcat('Histogram with Log-Normal Fit for the species ', {' '}, name_species));
% axis([0 20 0 inf]); % Adjust axis limits as needed
legend('Histogram', 'Log-Normal Fit');
end