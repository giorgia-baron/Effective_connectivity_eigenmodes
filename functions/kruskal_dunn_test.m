function kruskal_dunn_test(data)
    % Function to perform Kruskal-Wallis test and Dunn’s post hoc test
    % Input: data - A 143x3 matrix where each column represents a group

    % Validate input
    if size(data,2) ~= 3
        error('Input matrix must have 3 columns (representing 3 groups).');
    end

    % Perform Kruskal-Wallis test directly on matrix columns
    [p,~,stats] = kruskalwallis(data, [], 'off');

    % Display results
    fprintf('Kruskal-Wallis test p-value: %.5f\n', p);
    
    % Check if Kruskal-Wallis test is significant
    if p < 0.01
        disp('Kruskal-Wallis test is significant, performing post hoc test...');

        % Perform Dunn’s post hoc test with Bonferroni correction
        c = multcompare(stats, 'CType', 'bonferroni');

        % Display results in table format
        results_table = array2table(c, 'VariableNames', ...
            {'Group1', 'Group2', 'Lower_CI', 'Diff', 'Upper_CI', 'P_Value'});

        disp('Pairwise comparisons (Dunn’s test with Bonferroni correction):');
        disp(results_table);
    else
        disp('No significant difference between groups.');
    end
end
