function plot_ICC_vs_blockFreq(ICC_block_matrix, RSN_block_matrix_run1, RSN_block_matrix_run2)
% PLOT_ICC_VS_BLOCKFREQ Scatter plot of ICC block vs frequency block matrices
% with linear regression lines for each rr and two runs.
%
% Inputs:
%   ICC_block_matrix         - 3D matrix of ICC block matrices (blocks x blocks x rr)
%   RSN_block_matrix_run1    - 3D matrix of frequency block matrices for run 1
%   RSN_block_matrix_run2    - 3D matrix of frequency block matrices for run 2

% Define colors for each rr
color_range_3 = [0    0.4471    0.7412;   % blue
                 0.9294    0.6941    0.1255;   % orange
                 0.8510    0.3255    0.0980];  % red

num_rr = size(ICC_block_matrix,3);

%% Figure for Run 1
figure('Name','Run 1'); hold on;
for rr = 1:num_rr
    % Flatten matrices
    ICC_vec = ICC_block_matrix(:,:,rr);
    ICC_vec = ICC_vec(:);
    
    freq_vec = RSN_block_matrix_run1(:,:,rr);
    freq_vec = freq_vec(:);
    
    % Scatter plot
    scatter(freq_vec, ICC_vec, 36, color_range_3(rr,:), 'filled');
    
    % Linear fit
    p = polyfit(freq_vec, ICC_vec, 1);
    x_fit = linspace(min(freq_vec), max(freq_vec), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'Color', color_range_3(rr,:), 'LineWidth', 1.5);
end
xlabel('Frequency block'); ylabel('ICC block');
title('Run 1'); grid on; hold off;

%% Figure for Run 2
figure('Name','Run 2'); hold on;
for rr = 1:num_rr
    % Flatten matrices
    ICC_vec = ICC_block_matrix(:,:,rr);
    ICC_vec = ICC_vec(:);
    
    freq_vec = RSN_block_matrix_run2(:,:,rr);
    freq_vec = freq_vec(:);
    
    % Scatter plot
    scatter(freq_vec, ICC_vec, 36, color_range_3(rr,:), 'filled');
    
    % Linear fit
    p = polyfit(freq_vec, ICC_vec, 1);
    x_fit = linspace(min(freq_vec), max(freq_vec), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'Color', color_range_3(rr,:), 'LineWidth', 1.5);
end
xlabel('Frequency block'); ylabel('ICC block');
title('Run 2'); grid on; hold off;

end
