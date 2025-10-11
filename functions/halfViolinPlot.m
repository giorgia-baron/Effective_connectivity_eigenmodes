function halfViolinPlot(data1, data2, varargin)
% halfViolinPlot - Plots a left/right split violin plot comparing two distributions.
%
% Syntax: halfViolinPlot(data1, data2)
%         halfViolinPlot(data1, data2, 'ParameterName', ParameterValue, ...)
%
% Inputs:
%   data1 - Numeric vector for left-side distribution
%   data2 - Numeric vector for right-side distribution
%
% Optional Name-Value Pairs:
%   'Bandwidth' - Bandwidth for KDE (default: automatic)
%   'Colors'    - 1x2 cell array of colors (default: {'b','r'})
%
% Example:
%   x = randn(100,1);
%   y = randn(100,1) + 1;
%   halfViolinPlot(x, y, 'Colors', {[0.2 0.6 0.8], [0.8 0.2 0.4]});

% Parse inputs
p = inputParser;
addRequired(p, 'data1', @isnumeric);
addRequired(p, 'data2', @isnumeric);
addParameter(p, 'Bandwidth', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'Colors', {'b','r'}, @(x) iscell(x) && numel(x) == 2);
parse(p, data1, data2, varargin{:});
bw = p.Results.Bandwidth;
colors = p.Results.Colors;

% Kernel density estimation
[y1, x1] = ksdensity(data1, 'Bandwidth', bw);
[y2, x2] = ksdensity(data2, 'Bandwidth', bw);

% Normalize densities
y1 = y1 / max(y1);
y2 = y2 / max(y2);

% Mirror the violins
y1 = -y1; % left side

% Plot
figure; hold on;
fill(y1, x1, colors{1}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
fill(y2, x2, colors{2}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(zeros(1,2), [min([x1 x2]), max([x1 x2])], 'k--'); % vertical center line

xlabel('Density');
ylabel('Value');
title('Split Violin Plot Comparison');
legend({'True DAG','Null DAG'}, 'Location', 'best');
set(gca, 'YDir','normal');
hold off;
end
