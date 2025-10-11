clear all
clc
close all

%load data
load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/group_scatter.mat');
%including 0 Eim (matlab2022b)
load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/S_Sigma_HCP_110625.mat
addpath(genpath('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/FINAL_VERSION_CODE/functions/'))

N_rois=74; %subcorticals included
N_subs=length(S_Sigma_new_decomp);
N_nets=13;


% %get ROI colors
% group_col = [0.494117647058824   0.184313725490196   0.556862745098039;0.301960784313725   0.745098039215686   0.933333333333333;...
%     0.466666666666667   0.674509803921569   0.188235294117647;1     0     1;0.909803921568627   0.909803921568627   0.619607843137255;...
%     0.929411764705882   0.694117647058824   0.125490196078431;0.850980392156863   0.325490196078431   0.098039215686275;1 0 0.5; 0 0 0.5;...
%     1 0 0;0.75 0.75 0;0 1 0.5;1.0000    1.0000    0.0667];
% 
% color_rois=[];
% for rr=1:N_rois
%     color_rois(rr,:)=group_col(group_num(rr),:);
% end

%reorder labels and colors
new_order=[];
new_RSN_bounds(1)=1;
for nn=1:N_nets
    idx_net=find(group_num==nn);
    new_order=[new_order,idx_net];
    new_RSN_bounds(nn+1)=length(new_order)+1;
end



%color of each range
color_range=[0    0.4471    0.7412; 0.9294    0.6941    0.1255;0.8510    0.3255    0.0980];

%% visualize matrix S

% Number of struct elements
num_structs = numel(S_Sigma_new_decomp)-1;

% Preallocate the 4D array
S_modes_all = zeros(74, 74, 4, num_structs);

% Loop over each struct element
for ss = 1:num_structs
    S_modes_all(:,:,:,ss) = S_Sigma_new_decomp(ss).S_run1_norm;
end

% Remove the 2nd slice along the 3rd dimension
S_modes_all(:, :, 2, :) = [];

% Display size to confirm
disp('Size of S_modes_all:');
disp(size(S_modes_all));  % Should be [74 74 4 143]

S_modes_range=S_modes_all;


for rr = 1:3
    figure(rr)
    clf

    % Compute mean matrix
    mean_S(:,:,rr) = mean(squeeze(S_modes_range(new_order, new_order, rr, :)), 3);
    std_S(:,:,rr)  = std(squeeze(S_modes_range(new_order, new_order, rr, :)), [], 3);

    %% --- Matrix Plot ---
    subplot(2,1,1)
    imagesc(mean_S(:,:,rr))
    title(['S matrix range = ', num2str(rr)]);
    axis square;
    w_col = [1 1 1];
    col1 = color_range(rr,:);
    col2 = [0.5020 0.5020 0.5020];
    [grad1, ~] = colorGradient(w_col, col2, 200);
    [grad2, ~] = colorGradient(w_col, col1, 200);
    colormap([flipud(grad1); grad2])
    if rr == 1
        caxis([-0.005 0.005])
    else
        caxis([-0.01 0.01])
    end
    colorbar
    set(gca,'XTick',[],'YTick',[])
    hold on
    for jj = 1:length(new_RSN_bounds)
        xline(new_RSN_bounds(jj)-0.5, 'k')
        yline(new_RSN_bounds(jj)-0.5, 'k')
    end
    set(gca,'TickLength',[0 0])

    %% --- Bar Plot: Column Strength (sum) ---
    subplot(2,1,2)
    col_strength = sum(mean_S(:,:,rr), 1);  % sum (not abs)
    
    % Bar color matching current range color
    bar_vals = col_strength;
    bar_handle = bar(bar_vals, 'FaceColor', color_range(rr,:), 'EdgeColor', 'none');
    
    hold on
    for jj = 1:length(new_RSN_bounds)
        xline(new_RSN_bounds(jj)-0.5, 'k', 'LineWidth', 1);
    end
    hold off

    xlim([0.5 N_rois+0.5])
    ylim([-0.5 0.5])
    ylabel('Column Sum')
    set(gca,'XTick',[])
    box on

end



%% Compute RSN-level matrices (asymmetric i,i blocks)
N_nets_allsub=8;
RSN_block_matrix = zeros(N_nets_allsub, N_nets_allsub, 3);
group_num_allsub=group_num;
group_num_allsub(63:end)=8;

new_order_allsub=[];
for nn=1:N_nets_allsub
    idx_net=find(group_num_allsub==nn);
    new_order_allsub=[new_order_allsub,idx_net];
    
end

RSN_names={'VIS','SOM','DAN','SVAN','LIM','CON','DMN','SUB'};

for rr = 1:3
    RSN_block_matrix(:,:,rr) = computeRSNBlockMatrix(mean_S(:,:,rr), group_num_allsub, new_order_allsub, N_nets_allsub);
end

% Plot RSN-level block matrices
for rr = 1:3
    figure('Name', ['RSN Block Matrix Range ', num2str(rr)], 'Color', 'w');
    imagesc(RSN_block_matrix(:,:,rr));
    axis square;
    title(['RSN Block Matrix (Range ', num2str(rr), ')']);
    
    % Custom colormap
    col1 = color_range(rr,:);
    col2 = [0.5020 0.5020 0.5020];
    [grad1,~] = colorGradient([1 1 1], col2, 200);
    [grad2,~] = colorGradient([1 1 1], col1, 200);
    colormap([flipud(grad1); grad2])
     if rr==1
        caxis([-0.005 0.005])
    else
        caxis([-0.01 0.01])
    end
    colorbar;

    % Axes
    xticks(1:N_nets_allsub); yticks(1:N_nets_allsub);
    xlabel('RSN'); ylabel('RSN');
    set(gca,'TickLength',[0 0])
    
    % Tick labels
    xticks(1:N_nets_allsub); yticks(1:N_nets_allsub);
    xticklabels(RSN_names); yticklabels(RSN_names);
    xtickangle(45)
    set(gca,'TickLength',[0 0], 'FontSize', 10)
    
end

