clear all
clc
close all

load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/modes_class.mat


% Number of subjects and ranges
numSubjects = numel(Div_ranges);
numRanges = 4;  % ranges are labeled as 1, 2, 3, 4

% Initialize result matrices
mode_counts_run1 = zeros(numRanges, numSubjects);
mode_counts_run2 = zeros(numRanges, numSubjects);

for subj = 1:numSubjects
    % Get range labels for both runs
    run1_labels = Div_ranges(subj).div_range_run1;  % 1 × numModes
    run2_labels = Div_ranges(subj).div_range_run2;  % 1 × numModes

    % Count occurrences of each range (1–4) for run1
    for r = 1:numRanges
        mode_counts_run1(r, subj) = sum(run1_labels == r);
        mode_counts_run2(r, subj) = sum(run2_labels == r);
    end
end

% Save the result to MAT file
% save('mode_counts_divranges.mat', 'mode_counts_run1', 'mode_counts_run2');
