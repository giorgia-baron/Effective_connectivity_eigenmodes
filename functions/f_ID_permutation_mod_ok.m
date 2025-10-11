function [idiff_median_null, success_median_null, idiff_mean_null, success_mean_null] = f_ID_permutation_mod_ok(ID, n_subj, null_runs)

% Initialize outputs

success_null = zeros(1,null_runs);
idiff_null = zeros(1,null_runs);

success_median_null = zeros(1,null_runs);
idiff_median_null = zeros(1,null_runs);

success_mean_null = zeros(1,null_runs);
idiff_mean_null = zeros(1,null_runs);

for run = 1:null_runs
    % Permute both rows and columns
    ID_perm = ID(randperm(n_subj), randperm(n_subj));

    % Initialize vectors
    Iself = zeros(1, n_subj);
    Iothers = zeros(1, n_subj);
    Idiff_indiv = zeros(1, n_subj);
    succ_rate = zeros(1, n_subj);

    % Temporarily replace diagonals with NaN
    for s = 1:n_subj
        Iself(s) = ID_perm(s, s);  % self
        ID_perm(s, s) = nan;       % exclude self from mean calculation
    end

    % Compute Iothers and Idiff
    for s = 1:n_subj
        Iothers(s) = 0.5 * (nanmean(ID_perm(s, :)) + nanmean(ID_perm(:, s))');
        Idiff_indiv(s) = Iself(s) - Iothers(s);
    end

    % Compute success rate
    for s = 1:n_subj
        succ_rate(s) = sum(Iself(s) > ID_perm(s, :)) + sum(Iself(s) > ID_perm(:, s));
    end
    sr_indiv_perc = (succ_rate * 100 / (2 * n_subj - 2))';  % individual success %

    % Save metrics for this permutation run
 
    idiff_median_null(run) = median(Idiff_indiv);
    success_median_null(run) = median(sr_indiv_perc);

    idiff_mean_null(run) = mean(Idiff_indiv);
    success_mean_null(run) = mean(sr_indiv_perc);
end
end
