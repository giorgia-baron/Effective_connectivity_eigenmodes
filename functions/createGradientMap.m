function grad_range = createGradientMap(col1, nSteps)
% createGradientMap - Generates a colormap from white to gray to a given color.
%
% Syntax:
%   grad_range = createGradientMap(col1, nSteps)
%
% Inputs:
%   col1    - 1x3 RGB vector (e.g., [0.8, 0.2, 0.1])
%   nSteps  - Number of color steps in the gradient (default: 200)
%
% Output:
%   grad_range - nSteps x 3 colormap matrix

    if nargin < 2
        nSteps = 200;
    end

    % Define intermediate gray color
    w_col = [0.6510, 0.6510, 0.6510];

    % Step 1: Create gradient from white to gray
    n1 = round(nSteps / 2);
    white = [0.901960784313726   0.901960784313726   0.901960784313726];
    grad1 = [linspace(white(1), w_col(1), n1)', ...
             linspace(white(2), w_col(2), n1)', ...
             linspace(white(3), w_col(3), n1)'];

    % Step 2: Create gradient from gray to target color
    n2 = nSteps - n1;
    grad2 = [linspace(w_col(1), col1(1), n2)', ...
             linspace(w_col(2), col1(2), n2)', ...
             linspace(w_col(3), col1(3), n2)'];

    % Combine both gradients
    grad_range = [grad1; grad2];
end
