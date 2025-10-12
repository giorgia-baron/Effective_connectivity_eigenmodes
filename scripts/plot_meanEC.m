clear all
clc
close all

load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/DATA_EC_HCP.mat
load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/group_scatter.mat');

N_rois=74; %subcorticals included
N_subs=length(EC_MODES);
N_nets=13;


%get ROI colors
group_col = [0.494117647058824   0.184313725490196   0.556862745098039;0.301960784313725   0.745098039215686   0.933333333333333;...
    0.466666666666667   0.674509803921569   0.188235294117647;1     0     1;0.909803921568627   0.909803921568627   0.619607843137255;...
    0.929411764705882   0.694117647058824   0.125490196078431;0.850980392156863   0.325490196078431   0.098039215686275;1 0 0.5; 0 0 0.5;...
    1 0 0;0.75 0.75 0;0 1 0.5;1.0000    1.0000    0.0667];

group_num_basic=group_num;
group_num_basic(63:74)=8;

%reorder labels and colors
new_order=[];
new_RSN_bounds(1)=1;
for nn=1:N_nets
    idx_net=find(group_num_basic==nn);
    new_order=[new_order,idx_net];
    new_RSN_bounds(nn+1)=length(new_order)+1;
end


% Extract effective connectivity matrices
EC_all = cat(3, EC_MODES.EC1);  % EC1 is assumed to be a 74x74 matrix per subject

% Reorder if needed
EC_all = EC_all(new_order, new_order, :);

% Compute mean and std across subjects
mean_EC = mean(EC_all, 3);
std_EC = std(EC_all, [], 3);

% --- Plotting ---
figure(1)
clf

%% --- Matrix Plot ---
% subplot(2,1,1)
imagesc(mean_EC)
title('Mean Effective Connectivity (EC)');
axis square;

% Define color gradient
w_col = [1 1 1];  % white
col1 = [0.3922    0.8314    0.0745];  % dark red or any color of choice
col2 = [0 0 1];  % gray
[grad1, ~] = colorGradient(w_col, col2, 200);
[grad2, ~] = colorGradient(w_col, col1, 200);
colormap([flipud(grad1); grad2])

% Adjust color scale as needed for EC
caxis([-0.05 0.05])
colorbar
set(gca,'XTick',[],'YTick',[])
hold on
for jj = 1:length(new_RSN_bounds)
    xline(new_RSN_bounds(jj)-0.5, 'k')
    yline(new_RSN_bounds(jj)-0.5, 'k')
end
set(gca,'TickLength',[0 0])

%% --- Bar Plot: Column Strength (sum) ---
% subplot(2,1,2)
% col_strength = sum(mean_EC, 1);  % outgoing influence from each node
% 
% % Bar color
% bar_vals = col_strength;
% bar_handle = bar(bar_vals, 'FaceColor', col1, 'EdgeColor', 'none');
% 
% hold on
% for jj = 1:length(new_RSN_bounds)
%     xline(new_RSN_bounds(jj)-0.5, 'k', 'LineWidth', 1);
% end
% hold off
% 
% xlim([0 N_rois+1])
% ylabel('Column Sum')
% set(gca,'XTick',[])
% box on
