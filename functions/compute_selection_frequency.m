function freq_map = compute_selection_frequency(mask_all)
% Plots selection frequency maps with node-level incoming and outgoing frequency
% and adds RSN boundary lines on the heatmap
%
% mask_all: [numROIs x numROIs x numSubjects x numRanges] binary matrix
% roi_labels: optional, cell array of ROI labels
% new_RSN_bounds: optional, vector of ROI indices indicating RSN boundaries

[~, ~, ~, numRanges] = size(mask_all);

    for rr = 1:numRanges
        % New figure for each range
%         figure('Name', sprintf('Selection Frequency - Range %d', rr), 'Position', [100, 100, 1000, 700]);
    
        % Get binary mask for this range
        mask = mask_all(:, :, :, rr);  % [numROIs x numROIs x numSubjects]
    
        % Frequency map (mean across subjects)
        freq_map(:, :, rr) = mean(mask, 3);  % [numROIs x numROIs]
    
%         % Outgoing: mean across rows (i.e., column-wise average)
%         out_freq = mean(freq_map(:, :, rr), 1);  % [1 x numROIs]
%     
%         % Incoming: mean across columns (i.e., row-wise average)
%         in_freq = mean(freq_map(:, :, rr), 2);   % [numROIs x 1]
    
       
    end

end
