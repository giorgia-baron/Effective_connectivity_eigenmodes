function RSN_block_matrix = computeRSNBlockMatrix(full_mat, group_num, new_order, N_nets)
% Computes RSN-level block matrix with separate upper/lower triangle means for diagonal blocks
% Inputs:
% - full_mat: NxN matrix (e.g., mean_S(:,:,rr))
% - group_num: 1xN vector with RSN labels for each ROI
% - new_order: 1xN vector with reordered ROI indices
% - N_nets: number of RSNs
% Output:
% - RSN_block_matrix: N_nets x N_nets matrix

    RSN_labels = 1:N_nets;
    RSN_block_matrix = nan(N_nets, N_nets);

    ordered_group_num = group_num(new_order);  % reorder group labels to match matrix

    for i = 1:N_nets
        idx_i = find(ordered_group_num == RSN_labels(i));
        for j = 1:N_nets
            idx_j = find(ordered_group_num == RSN_labels(j));
            block = full_mat(idx_i, idx_j);

%             if i ~= j
                RSN_block_matrix(i,j) = mean(block(:));  % between-RSN average
%             else
%                 % Within-RSN: split upper and lower
%                 upper_vals = block(triu(true(length(idx_i)), 1));
%                 lower_vals = block(tril(true(length(idx_i)), -1));
%                 RSN_block_matrix(i,j) = mean(upper_vals);
%                 RSN_block_matrix(j,i) = mean(lower_vals);
%             end
        end
    end
end
