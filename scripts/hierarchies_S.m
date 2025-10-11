


clear all
clc
close all

%% load data and paths

load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/group_scatter.mat');
load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/S_Sigma_HCP_110625.mat
load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/node_labels_Schaefer.mat
addpath(genpath('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/FINAL_VERSION_CODE/functions/'))


%get ROI colors
group_col = [0.494117647058824   0.184313725490196   0.556862745098039;0.301960784313725   0.745098039215686   0.933333333333333;...
    0.466666666666667   0.674509803921569   0.188235294117647;1     0     1;0.909803921568627   0.909803921568627   0.619607843137255;...
    0.929411764705882   0.694117647058824   0.125490196078431;0.850980392156863   0.325490196078431   0.098039215686275;1 0 0.5; 0 0 0.5;...
    1 0 0;0.75 0.75 0;0 1 0.5;1.0000    1.0000    0.0667];
group_net=unique(group,'stable');

N_rois=74; %subcorticals included
N_subs=143;
N_nets=13;


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

%ROI colors
color_rois=[];
for rr=1:N_rois
    color_rois(rr,:)=group_col(group_num(rr),:);
end


% load DAGs and associated S matrices for each range
dag_r1=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_run1_range1.mat').DAG_all_range1;
dag_r2=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_run1_range2.mat').DAG_all_range2;
dag_r3=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_run1_range3.mat').DAG_all_range3;

Sr1=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_run1_range1.mat').MATRIX_S_range1;
Sr2=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_run1_range2.mat').MATRIX_S_range2;
Sr3=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_run1_range3.mat').MATRIX_S_range3;

N_subs=size(Sr1,3);

%% extract DAG for each subject and range

for rr_dag=1:3

    count_subj=1;
    DAG_sparse_dg_all=struct();
    
    for ss=1:N_subs
       
        %create DAG mask to retain only stable links 
        DAG_bin=not(eval(['dag_r',num2str(rr_dag)]).DAG(:,:,:,ss)==0); %selezionare range
        DAG_bin_sum=sum(DAG_bin,3);
   
    
        %keep only links present in all iterations
        mask_DAG=double(DAG_bin_sum==100);
        

        MATRIX_S=eval(['Sr',num2str(rr_dag)]);
        DAG_sparse(:,:,ss)=MATRIX_S(:,:,ss).*mask_DAG; 
        DAG_frequency_mat(:,:,ss,rr_dag)=mask_DAG;

        if sum(DAG_sparse(:,:,ss)<0,'all')
            disp('neg values detected')
            pause
        end

        DAG_sparse_dg=digraph(DAG_sparse(:,:,ss),node_labels);
        if isdag(DAG_sparse_dg) 
            DAG_sparse_dg_all(count_subj).DAG_sparse=DAG_sparse_dg;
            count_subj=count_subj+1;
        else
            disp('detected no DAG graph')
            pause
        end
    
        % get the full DAG (with both positive and negative values) and
        % store for plot 
        mask_DAG_tmp=mask_DAG;
        % Extract lower triangular part (excluding diagonal)
        L = tril(mask_DAG_tmp, -1);   
        % Transpose the lower triangular part
        L_T = L';
        % Extract the upper triangular part 
        U = triu(mask_DAG_tmp,1);
        
        % Sum the transposed lower part with the upper part
        mask_DAG_fin = U + L_T;
        mask_DAG_fin=mask_DAG_fin+mask_DAG_fin';
        
    %     figure,imagesc(mask_DAG_fin)
        mask_all_ord_freq_selection(:,:,ss,rr_dag)=mask_DAG_fin(new_order,new_order);

    end
    

    
   
    %% topological sorting of DAGs

    for ss=1:N_subs
        
    
        clear rank_toposort
        DAG_sparse_sub=DAG_sparse_dg_all(ss).DAG_sparse;
        outD(ss,:)=outdegree(DAG_sparse_sub); %controllare se i link con degree zero sono alla fine
        zero_outD(ss).idx_zero_out=find(outD(ss,:)==0);%regioni con outdeg=0
        inD(ss,:)=indegree(DAG_sparse_sub);
        zero_inD(ss).idx_zero_in=find(inD(ss,:)==0); %regioni con indeg=0
        zero_D(ss).idx_zero_all=intersect(zero_outD(ss).idx_zero_out, zero_inD(ss).idx_zero_in); %regioni che hanno un degree totale =0 -> mettere a NaN
         if sum(zero_D(ss).idx_zero_all)>0
            disp('zero deg nodes detected')
            pause
         end

        %to order according to the outdegree (problem of applying coherently toposort)-> mi ritorna le regioni da quella con outdegree più alto a più basso (alla fine mette quelle con outdeg=0)
        [~, idx_outdegree_old] = sort(DAG_sparse_sub.outdegree, 1, "descend");

        outdeg = DAG_sparse_sub.outdegree;   % out-degree
        outstr = sum(DAG_sparse(:,:,ss),2); % out-strength (sum of outgoing weights)
        % Create a 2-column matrix: [out-degree, out-strength]
        sortMat = [outdeg(:), outstr(:)];
        
        % Sort by first column descending, then second column descending
        [~, idx_outdegree] = sortrows(-sortMat, [1 2]);  


        idx_sort_outdgree_sorted = toposort(reordernodes(DAG_sparse_sub, idx_outdegree), "Order", "stable"); %ordina in base ai nodi che sono più source verso quelli più sink
        idx_toposort = idx_outdegree(idx_sort_outdgree_sorted); %questo mi serve per vedere, per ogni soggetto, l'ordine topologico delle ROI (non ordinate per network) 
        [~,rank_toposort]=sort(idx_toposort); %questo mi serve per ottenere gli indici dell'ordine topologico considerando l'ordine delle network nell'atlente (da VIS a DMN e poi SUBCORTICAL) - quindi mi dirà che il VIS (prime due posizioni sono nell'ordine topologico 54-55, quindi sono più verso l'essere sink piuttosto che source) 
         
        if not(isequal(sort(idx_toposort)', 1:N_rois))
            disp('Toposort did not return a full permutation of ROIs')
            pause
        end

        rank_toposorted_nozeroDeg=rank_toposort;
        rank_toposorted_nozeroDeg(zero_D(ss).idx_zero_all)=NaN;
        
        %save ranking for each subject
        asym_rank_all_nozeroDeg(ss,:)=rank_toposorted_nozeroDeg;
   
    end
    
    
    %% Heatmap sorted by descendent ranking
    
    [heatmap_rank_nonzero_ord,idx_subj_order]=sort(asym_rank_all_nozeroDeg,'descend');
    [num_subjects, num_columns] = size(heatmap_rank_nonzero_ord);

    sub_idx=1:N_subs;
    sub_idx_cell=num2cell(sub_idx);
    xLabels = node_labels(new_order); 
    
    hh=figure
    heat=heatmap(heatmap_rank_nonzero_ord(:,new_order))
    colormap jet
    colorbar
    xlabel('ROIs')
    ylabel('Subjects')
    heat.YDisplayLabels = repmat({''}, size(heatmap_rank_nonzero_ord(:,new_order), 1), 1); % Rimuove le etichette Y
    heat.XDisplayLabels = xLabels; 
    heat.Title='Single subject hierachical rankings';
    heat.FontSize=7;
    set(hh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    %save for ball and sticks plot
    ascend_box_ord(:,rr_dag)=mean(asym_rank_all_nozeroDeg);

    %% RSN ranking in violin plots

    asym_rank_all_nozeroDeg_vec=asym_rank_all_nozeroDeg';
    asym_rank_all_nozeroDeg_vec=asym_rank_all_nozeroDeg_vec(:);
    group_rep=repmat(group,1,N_subs);
    group_num_rep=repmat(group_num,1,N_subs);
    rsn_rank=grpstats(asym_rank_all_nozeroDeg_vec,group_rep',@(x)mean(x,1));
    [~,rsn_rank_ord_idx]=sort(rsn_rank,'ascend');
    
    %save for statistical trend test
    jtrend_ranking(:,:,rr_dag)=asym_rank_all_nozeroDeg;
    [~,rsn_rank_ranges(:,rr_dag)]=sort(rsn_rank_ord_idx,'ascend');
    rsn_rank_mean(:,rr_dag)=rsn_rank;

    %order ROIs according to RSN ascending order
    ascend_rois_rsn_order=[];
    %jtrend_label=zeros(N_subs,N_rois);-->ERROR GB
    for ll=1:N_nets
    
        idx_net_box=group_num==rsn_rank_ord_idx(ll);
        jtrend_label_rr(:,idx_net_box)=rsn_rank_ranges(ll,rr_dag);
        rank_within_curr_net=mean(asym_rank_all_nozeroDeg(:,idx_net_box));
        [~,ord_within_curr_net_idx]=sort(rank_within_curr_net,'ascend');
        idx_net_sort=find(idx_net_box);
        idx_net_sort=idx_net_sort(ord_within_curr_net_idx);
        ascend_rois_rsn_order=[ascend_rois_rsn_order, idx_net_sort];
        
    end
    
    jtrend_label(:,:,rr_dag)=repmat(jtrend_label_rr,143,1);

    %violin plot
    color_rois_ascend_box=color_rois(ascend_rois_rsn_order,:);
    f=figure
    vplot=violinplot(asym_rank_all_nozeroDeg(:,ascend_rois_rsn_order),node_labels(ascend_rois_rsn_order),...
    color_rois_ascend_box,'ViolinColor',[0 0 0],...
    'ShowData',false,'GroupOrder',group(ascend_rois_rsn_order),'Width',0.4,'BoxColor',[0 0 0],...
    'MedianColor',[0 0 0],'ViolinAlpha',0.6,'ShowMean',true)
    ylim([1 N_rois])
    axis tight;
    grid on;
    hold off;
    ylabel('Hierarchical ranking')
    title(['S graph R_', num2str(rr_dag)])
    subtitle('ROI level ranking','FontSize',8)
    set(f, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(gca, 'Position', [0.05, 0.1, 0.9, 0.8]); % Expand axes inside figure
   
    
%     pause
    close all

end

%% plot DAG frequency selection at RSN level

freq_map=compute_selection_frequency(mask_all_ord_freq_selection);

N_nets_allsubcort=8;
group_num_allsubcort=group_num;
group_num_allsubcort(63:end)=8;

%block matrix
RSN_block_matrix = zeros(N_nets_allsubcort, N_nets_allsubcort, 3);
RSN_names={'VIS','SOM','DAN','SVAN','LIM','CON','DMN','SUB'};

new_order_allsubcort=[];
for nn=1:N_nets
    idx_net=find(group_num_allsubcort==nn);
    new_order_allsubcort=[new_order_allsubcort,idx_net];
    
end


% Plot RSN-level block matrices
for rr = 1:3

    RSN_block_matrix(:,:,rr) = computeRSNBlockMatrix(freq_map(:,:,rr), group_num_allsubcort, new_order_allsubcort, N_nets_allsubcort);

    figure('Name', ['RSN DAG Block Matrix Range ', num2str(rr)], 'Color', 'w');
    imagesc(RSN_block_matrix(:,:,rr));
    axis square;
    title(['RSN DAG Block Matrix (Range ', num2str(rr), ')']);
    
    % Custom colormap
    col1 = color_range(rr,:);
    col2 = [0.5020 0.5020 0.5020];
    [grad1,~] = colorGradient([1 1 1], col2, 200);
    [grad2,~] = colorGradient([1 1 1], col1, 200);
    colormap([flipud(grad1); grad2])
    colorbar;

    % Axes
    xticks(1:N_nets_allsubcort); yticks(1:N_nets_allsubcort);
    xlabel('RSN'); ylabel('RSN');
    set(gca,'TickLength',[0 0])
    
    % Tick labels
    xticks(1:N_nets_allsubcort); yticks(1:N_nets_allsubcort);
    xticklabels(RSN_names); yticklabels(RSN_names);
    xtickangle(45)
    set(gca,'TickLength',[0 0], 'FontSize', 10)
    
end


%% network spider plot across ranges

P = rsn_rank_ranges';

% Delete variable in workspace if exists
if exist('s', 'var')
    delete(spid_plot);
end

% Spider plot
figure
spid_plot = spider_plot_class(P);
spid_plot.AxesWebType = 'circular';
spid_plot.FillOption = 'on';
spid_plot.AxesLabels= group_net;
spid_plot.LegendLabels={'R1','R2','R3'};
spid_plot.Color=color_range;
max_p=max(P,[],'all');
spid_plot.AxesInterval = 4;
spid_plot.AxesPrecision = 2;
spid_plot.AxesDisplay = 'one';
spid_plot.AxesLabelsEdge = 'none'; 
spid_plot.AxesStart= 1.74;
spid_plot.MarkerSize=45;
spid_plot.Direction='counterclockwise';
spid_plot.AxesLimits = [0, 0, 0, 0, 0,0,0,0,0,0,0,0,0; max_p,max_p,max_p,max_p,max_p,max_p,max_p,max_p,max_p,max_p,max_p,max_p,max_p];

%% save group-level DAG for brain rendering

for rr_dag=1:3
    
    N_DAGs=size(DAG_frequency_mat,3);
    DAG_frequency_mat_perc=sum(DAG_frequency_mat(:,:,:,rr_dag),3)./N_DAGs*100;

    %extract edge matrix
    prc_99=prctile(DAG_frequency_mat_perc(:),99);
    DAG_edge(:,:,rr_dag)=DAG_frequency_mat_perc.*(DAG_frequency_mat_perc>=prc_99);
    
end

%% Trend test on ranking

for rr=1:3

    data=jtrend_ranking(:,:,rr);
    g=jtrend_label(:,:,rr);
    jttrend([data(:),g(:)])
    
end

%% compare number of removed cycles null vs true model

figure;  % Create a single figure for all subplots

for rr = 1:3
    % Load data
    load(['/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/Results_test_null_DAG_run1_range', num2str(rr), '.mat'])

    % ----------- SUM OF WEIGHTS OF REMOVED EDGES -------------
    WR_tot = eval(['WR_range', num2str(rr)]);
    WR_rand_tot = eval(['WR_rand_range', num2str(rr)]);
    
    if sum(isnan(WR_rand_tot),'all')>0 
        disp('NaN values detected in random WR'), pause
        rr
    end
    if sum(isnan(WR_tot),'all')>0 %modif GB
        sum(isnan(WR_tot))
        WR_tot(isnan(WR_tot))=0;
    end
    
    % subplot(3, 2, 2 * rr - 1);  % Left column (odd indices)
    halfViolinPlot(mean(WR_tot)', mean(WR_rand_tot)')
    title(['Sum Weights Removed - Range ', num2str(rr)])
    ylabel('Sum of Weights')
    xticklabels({'Original', 'Randomized'})

    % ----------- NUMBER OF REMOVED EDGES -------------
    ER_tot = eval(['ER_range', num2str(rr)]);
    ER_rand_tot = eval(['ER_rand_range', num2str(rr)]);
    
    if sum(isnan(ER_rand_tot),'all')>0 
        disp('NaN values detected in random ER'), pause
        rr
    end
    if sum(isnan(ER_tot),'all')>0
        sum(isnan(ER_tot))
        ER_tot(isnan(ER_tot))=0;
    end

    % subplot(3, 2, 2 * rr);  % Right column (even indices)
    halfViolinPlot(mean(ER_tot)', mean(ER_rand_tot)');
    title(['Number of Edges Removed - Range ', num2str(rr)])
    ylabel('Number of Edges')
    xticklabels({'Original', 'Randomized'})

    % ----------- STATISTICAL TESTING -------------
     data1e = mean(ER_tot)';
    data2e = mean(ER_rand_tot)';
    data1w = mean(WR_tot)';
    data2w = mean(WR_rand_tot)';
    
    % [h1, p1] = lillietest(data1);
    % [h2, p2] = lillietest(data2);
    % [h_var, p_var] = vartest2(data1, data2);
    [h_e, p_e] = ranksum(data1e, data2e)  % Equal variance
    [h_w, p_w] = ranksum(data1w, data2w)  % Equal variance
    % [h_uneq, p_uneq] = ttest2(data1, data2, 'Vartype','unequal');  % Welch

    pause
end

% sgtitle('Run 1: Removed Edge Analysis (Ranges 1-3)')

%% ball and sticks plot for range hierarchy

%load atlas nodes coordinates
load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Schaefer_100Parcels_7Networks_AAL3_2_MNI152FSL_2mm_clustered_coords_matrix.mat

coord=a;
%switch cerebellum and hippocampus in LH
coord=coord([1:66,68,67,69:74],:);

saveFolder='/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/';

for rr_dag=1:3
    
    fileNameNode=['DAG_run1_node' num2str(rr_dag) '.node'];
    fileNameEdge=['DAG_run1_edge' num2str(rr_dag) '.edge'];
    tmp=coord;
    %RSN class
    tmp(:,4)=group_num;
    %ball size proportional to group-mean ranking
    tmp(:,5)=1./ascend_box_ord(:,rr_dag);
    dlmwrite([saveFolder,fileNameNode],tmp,'\t');

    %from graphic visualization seems DAG to be transposed
    dlmwrite([saveFolder,fileNameEdge],DAG_edge(:,:,rr_dag)','\t');

end

addpath(genpath('/nfsd/biopetmri4/Users/Giorgia/Softwares/NITRC-multi-file-downloads(1)'))

