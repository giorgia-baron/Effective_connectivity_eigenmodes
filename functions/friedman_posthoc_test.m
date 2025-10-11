function friedman_posthoc_test(data,pval_thr)
    % Function to perform Friedman test and post hoc multiple comparisons
    % Input: data - An N x 3 matrix where each row is a subject,
    % and each column is a condition (e.g., R1, R2, R3)

    % Validate input
    if size(data,2) ~= 3
        error('Input matrix must have 3 columns (representing 3 conditions).');
    end

    % Number of subjects
    nSubjects = size(data,1);

    % Group labels for each condition
    group_labels = [1 2 3]; % Corresponds to R1, R2, R3

    % Friedman test (repeated measures)
    [p, tbl, stats] = friedman(data, 1, 'off');

    fprintf('Friedman test p-value: %.5f\n', p);

    if p < pval_thr
        disp('Friedman test is significant, performing post hoc comparisons...');

        % Post hoc multiple comparisons with Bonferroni correction
        posthoc = multcompare(stats, 'CType', 'bonferroni');

        % Create readable table
        results_table = array2table(posthoc, 'VariableNames', ...
            {'Group1', 'Group2', 'Lower_CI', 'Mean_Diff', 'Upper_CI', 'P_Value'});

        % Label the groups as R1, R2, R3
        group_names = {'R1', 'R2', 'R3'};
        results_table.Group1 = group_names(results_table.Group1)';
        results_table.Group2 = group_names(results_table.Group2)';

        disp('Pairwise comparisons (Bonferroni-corrected):');
        disp(results_table);
    else
        disp('No significant difference between conditions.');
    end
end
