function plot_vertical_distributions(dist1, dist2, edges, varargin)
% PLOT_VERTICAL_DISTRIBUTIONS Plots two distributions attached vertically.
%   plot_vertical_distributions(dist1, dist2, edges)
%   dist1, dist2: data vectors or histograms
%   edges: bin edges to use (shared for both)
%   Optional Name-Value pairs:
%       'Normalize' (true/false, default: true)
%       'Labels' (cell array with two labels)
%       'Colors' (cell array with two color specs)

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'Normalize', true, @islogical);
    addParameter(p, 'Labels', {'Distribution 1', 'Distribution 2'});
    addParameter(p, 'Colors', {'b', 'r'});
    parse(p, varargin{:});
    normalize = p.Results.Normalize;
    labels = p.Results.Labels;
    colors = p.Results.Colors;

    % Compute histograms if raw data is provided
    if numel(dist1) ~= numel(edges)-1
        dist1 = histcounts(dist1, edges);
    end
    if numel(dist2) ~= numel(edges)-1
        dist2 = histcounts(dist2, edges);
    end

    % Normalize if required
    if normalize
        dist1 = dist1 / max(dist1);
        dist2 = dist2 / max(dist2);
    end

    % Bin centers
    centers = edges(1:end-1) + diff(edges)/2;

    % Plot
    figure; hold on;
    bar(centers, dist1, 'FaceColor', colors{1}, 'EdgeColor', 'none');
    bar(centers, -dist2, 'FaceColor', colors{2}, 'EdgeColor', 'none');

    % Formatting
    yline(0, '--k');
    xlabel('Value');
    ylabel('Normalized Count');
    legend(labels);
    title('Vertical Comparison of Two Distributions');
    grid on;

    % Adjust Y limits
    yMax = max([dist1(:); dist2(:)]);
    ylim([-yMax, yMax] * 1.1);
end
