function raincloud_plot(data)
% Remove NaN values from data
data = data(~any(isnan(data), 2), :);  % Remove rows with any NaN values

% Number of distributions (columns in the data matrix)
num_vars = size(data, 2);

% Set up the figure for the plot
figure;
hold on;

% Define colors for each distribution
colors = [0    0.4471    0.7412; 0.9294    0.6941    0.1255;0.8510    0.3255    0.0980];  

% Overlay scatter points and smoothed density plots for each column
for i = 1:num_vars
    % Get the current variable data (without NaNs)
    current_data = data(:, i);
    
    % Create random jitter for x-axis position to spread the scatter points
    jitter = (rand(size(current_data)) - 0.5) * 0.2 + i; % Small random offset around i
    
    % Plot scatter points
    scatter(jitter, current_data, 20, 'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.5);

    % Compute histogram for density estimation
    [f, edges] = histcounts(current_data, 'Normalization', 'pdf', 'BinWidth', 0.05);
    
    % Get bin centers from histogram edges
    xi = edges(1:end-1) + diff(edges) / 2;
    
    % Smooth the density using interpolation
    xi_smooth = linspace(min(xi), max(xi), 200);  % Interpolate over a finer range
    f_smooth = interp1(xi, f, xi_smooth, 'pchip', 0);  % Smooth using piecewise cubic interpolation
    
    % Normalize density plot for proper scaling
    f_smooth = f_smooth / max(f_smooth) * 0.3;  % Scale width of density plots

    % Construct filled density plot
    fill(i + [0 f_smooth 0], [xi_smooth(1) xi_smooth xi_smooth(end)], colors(i, :), ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none');
end

% Adjust axis properties
xlim([0, num_vars + 1.5]);  % Extend x-axis to fit density plots
ylim([min(data(:)), max(data(:))]);
set(gca, 'XTick', 1:num_vars, 'XTickLabel', {'r1', 'r2', 'r3'});
xlabel('Variables');
title('Raincloud Plot: Smoothed Distribution with Scatter Overlay');
hold off;
end
