function plot_ICC_vs_freq(freq_map_run1, freq_map_run2, ICC_freq, triu_mask)
% PLOT_ICC_VS_FREQ Plot ICC boxplots with scatter points based on frequency percentiles
%
% Inputs:
%   freq_map_run1 - 3D matrix of frequency maps for run 1 (channels x channels x rr)
%   freq_map_run2 - 3D matrix of frequency maps for run 2
%   ICC_freq      - 3D matrix of ICC values (channels x channels x rr)
%   triu_mask     - logical mask for upper triangular part of matrices

num_rr = size(freq_map_run1, 3);  % number of rr

%% Figure for Run 1
figure('Name','Run 1'); hold on;
for rr = 1:num_rr
    % Upper triangular values
    freq1_triu = freq_map_run1(:,:,rr);
    freq1_triu = freq1_triu(triu_mask);
    
    ICC_triu = ICC_freq(:,:,rr);
    ICC_triu = ICC_triu(triu_mask);
    
    % Boxplot
    boxplot(ICC_triu,'Positions',rr,'Colors','k','Widths',0.5);
    
    % 50th percentile threshold
    p90 = prctile(freq1_triu,99);
    p10= prctile(freq1_triu,1);
    
    % Scatter points
    x_jitter = rr + 0.1*randn(size(ICC_triu));
    scatter(x_jitter(freq1_triu <= p10), ICC_triu(freq1_triu <= p10), 36, 'r','filled');
    scatter(x_jitter(freq1_triu >= p90), ICC_triu(freq1_triu >= p90), 36, 'b','filled');
end
xlabel('rr'); ylabel('ICC'); title('Run 1'); hold off;

%% Figure for Run 2
figure('Name','Run 2'); hold on;
for rr = 1:num_rr
    freq2_triu = freq_map_run2(:,:,rr);
    freq2_triu = freq2_triu(triu_mask);
    
    ICC_triu = ICC_freq(:,:,rr);
    ICC_triu = ICC_triu(triu_mask);
    
    % Boxplot
    boxplot(ICC_triu,'Positions',rr,'Colors','k','Widths',0.5);
    
    % Use run1 percentile threshold for coloring
    freq2_triu = freq_map_run2(:,:,rr);
    freq2_triu = freq2_triu(triu_mask);
    p90 = prctile(freq2_triu,99);
    p10= prctile(freq2_triu,1);
    
    % Scatter points
    x_jitter = rr + 0.1*randn(size(ICC_triu));
    scatter(x_jitter(freq2_triu <= p10), ICC_triu(freq2_triu <= p10), 36, 'r','filled');
    scatter(x_jitter(freq2_triu >= p90), ICC_triu(freq2_triu >= p90), 36, 'b','filled');
end
xlabel('rr'); ylabel('ICC'); title('Run 2'); hold off;

end
