function densityScatterChart_new(x, y, numBins, cmap,cbarLimits)
    % densityScatterChart: Create a density scatter plot where color
    % represents the density of points in a 2D space.
    %
    % INPUTS:
    %   x - A vector of x-coordinates of the points.
    %   y - A vector of y-coordinates of the points.
    %   numBins - (optional) The number of bins to use for the 2D histogram. Default is 50.
    %   cmap - (optional) The colormap to use. If not provided, default is 'jet'.

    if nargin < 3
        numBins = 50;  % Default number of bins if not specified
    end
    if nargin < 4
        cmap = 'jet';  % Default colormap if not specified
    end

    % Step 1: Calculate the 2D histogram of the points to get density values
    [N, Xedges, Yedges] = histcounts2(x, y, numBins);
    
    % Step 2: Compute the bin centers (midpoints of bins)
    Xcenters = (Xedges(1:end-1) + Xedges(2:end)) / 2;
    Ycenters = (Yedges(1:end-1) + Yedges(2:end)) / 2;
    
    % Step 3: Assign each point to its respective bin in the 2D histogram
    [~, xBin] = histc(x, Xedges);
    [~, yBin] = histc(y, Yedges);
    
    % Step 4: Get the density of each point based on its bin
    density = arrayfun(@(i) N(xBin(i), yBin(i)), 1:length(x));
    
    % Step 5: Plot the scatter plot with colors corresponding to density
    figure;
    scatter(x, y, 20, density, 'filled');  % 20 is the marker size
    colorbar;  % Show a colorbar to indicate the density scale
    colormap(cmap);  % Set the colormap
    caxis(cbarLimits);  % Set colorbar limits
    title('Density Scatter Plot');
    xlabel('X');
    ylabel('Y');
end