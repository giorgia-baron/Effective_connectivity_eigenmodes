clear all
clc
close all


%% load useful data

addpath(genpath('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/FINAL_VERSION_CODE/functions'))
addpath(genpath('/nfsd/biopetmri4/Users/Giorgia/Softwares/NIfTI_20140122/'))

load '/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/group_scatter.mat'
load '/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/node_labels_Schaefer.mat'

N_rois=74;
N_nets=13;
N_subs=143;

%build upper triangular mask
mask_ut = logical(triu(ones(N_rois,N_rois),1));
triu_mask=find(mask_ut);

%get ROI colors
group_col = [0.494117647058824   0.184313725490196   0.556862745098039;0.301960784313725   0.745098039215686   0.933333333333333;...
    0.466666666666667   0.674509803921569   0.188235294117647;1     0     1;0.909803921568627   0.909803921568627   0.619607843137255;...
    0.929411764705882   0.694117647058824   0.125490196078431;0.850980392156863   0.325490196078431   0.098039215686275;1 0 0.5; 0 0 0.5;...
    1 0 0;0.75 0.75 0;0 1 0.5;1.0000    1.0000    0.0667];

%define ICC threshold as in ref: https://www.nature.com/articles/s41598-018-25089-1
icc_mat_thr=80; 


%reorder labels and colors
new_order=[];
new_RSN_bounds(1)=1;
for nn=1:N_nets
    idx_net=find(group_num==nn);
    new_order=[new_order,idx_net];
    new_RSN_bounds(nn+1)=length(new_order)+1;
end

%set color of each range
color_range=[0    0.4471    0.7412; 0.9294    0.6941    0.1255; 0.8510    0.3255    0.0980];% blu, arancione, rosso 

%% load DAGs for both runs

dag_run1_r1=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_RUN1_range1.mat').DAG_all_range1;
dag_run1_r2=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_RUN1_range2.mat').DAG_all_range2;
dag_run1_r3=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_RUN1_range3.mat').DAG_all_range3;

dag_run2_r1=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN2/DAG_data_RUN2_range1.mat').DAG_all_range1;
dag_run2_r2=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN2/DAG_data_RUN2_range2.mat').DAG_all_range2;
dag_run2_r3=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN2/DAG_data_RUN2_range3.mat').DAG_all_range3;


%% create TEST and RETEST mat for DAGs

path_S='/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA';
load(fullfile(path_S,'S_Sigma_HCP_110625.mat'))

count=1;
for subj=1:N_subs
    
   
    if not(isempty(S_Sigma_new_decomp(subj).S_run1)) && not(isempty(S_Sigma_new_decomp(subj).S_run2))
        
        %RUN1
        S_sub_all_run1=S_Sigma_new_decomp(subj).S_run1;
        if size(S_sub_all_run1,3)==4
            S_sub_all_run1(:,:,2)=[];
            if (S_Sigma_new_decomp(subj).S_run1(:,:,1)==S_sub_all_run1(:,:,1)) & (S_Sigma_new_decomp(subj).S_run1(:,:,3)==S_sub_all_run1(:,:,2)) & (S_Sigma_new_decomp(subj).S_run1(:,:,4)==S_sub_all_run1(:,:,3))
                disp('Range 2 removed correctly!')
            end
        end
        
        %RUN2
        S_sub_all_run2=S_Sigma_new_decomp(subj).S_run2;
        if size(S_sub_all_run2,3)==4
        S_sub_all_run2(:,:,2)=[];
    
            if (S_Sigma_new_decomp(subj).S_run2(:,:,1)==S_sub_all_run2(:,:,1)) & (S_Sigma_new_decomp(subj).S_run2(:,:,3)==S_sub_all_run2(:,:,2)) & (S_Sigma_new_decomp(subj).S_run2(:,:,4)==S_sub_all_run2(:,:,3))
                disp('Range 2 removed correctly!')
            end
        end

        test_subj=S_sub_all_run1;
        retest_subj=S_sub_all_run2;
    else
        continue;
    end

    %check both runs are present
    if not(isempty(test_subj)) && not(isempty(retest_subj))
       
        %split ranges
        for rr=1:3
            
            %RUN1
            tmp_test=test_subj(:,:,rr);

            DAG_bin=not(eval(['dag_run1_r',num2str(rr)]).DAG(:,:,:,subj)==0); %selezionare range
            DAG_bin_sum=sum(DAG_bin,3);
            mask_DAG=DAG_bin_sum==100;

            %get corresponding mask for upper triangular part only

            %transpose matrix so source=columns
            mask_DAG=mask_DAG';
            % Extract lower triangular part (excluding diagonal)
            L = tril(mask_DAG, -1); 
            % Transpose the lower triangular part
            L_T = L';  
            % Extract the upper triangular part 
            U = triu(mask_DAG,1);
            
            % Sum the transposed lower part with the upper part
            mask_DAG_fin = U + L_T;
            mask_DAG_fin_run1=mask_DAG_fin+mask_DAG_fin';

            %sanity checks
%             figure,subplot(121),imagesc(mask_DAG_fin),subplot(122),imagesc(mask_DAG),pause
            if length(unique(mask_DAG_fin_run1))>2
                disp('error in building triangular mask')
                pause
            end
            mask_all_ord_freq_selection_run1(:,:,subj,rr)=mask_DAG_fin_run1(new_order,new_order);


            tmp_test=tmp_test-diag(diag(tmp_test));
            tmpt_triu=tmp_test(triu_mask);

            TEST(:,count,rr)=tmpt_triu(:);


            %RUN2
            tmp_retest=retest_subj(:,:,rr);
  
            DAG_bin=not(eval(['dag_run2_r',num2str(rr)]).DAG(:,:,:,subj)==0); %selezionare range
            DAG_bin_sum=sum(DAG_bin,3);
            mask_DAG=DAG_bin_sum==100;
            
            %get corresponding mask for upper triangular part only

            %transpose matrix so source=columns
            mask_DAG=mask_DAG';
            % Extract lower triangular part (excluding diagonal)
            L = tril(mask_DAG, -1); 
            % Transpose the lower triangular part
            L_T = L';  
            % Extract the upper triangular part 
            U = triu(mask_DAG,1);
            
            % Sum the transposed lower part with the upper part
            mask_DAG_fin = U + L_T;
            mask_DAG_fin_run2=mask_DAG_fin+mask_DAG_fin';

            %sanity checks
%             figure,imagesc(mask_DAG_fin),pause
            if length(unique(mask_DAG_fin_run2))>2
                disp('error')
                pause
            end

            mask_all_ord_freq_selection_run2(:,:,subj,rr)=mask_DAG_fin_run2(new_order,new_order);


            tmp_retest=tmp_retest-diag(diag(tmp_retest));
            tmpr_triu=tmp_retest(triu_mask);
            RETEST(:,count,rr)=tmpr_triu(:);
        end

        count=count+1;

    end

end

%check that count matches N_subs
count==N_subs+1

%compute selection frequency matrices
freq_map_run1=compute_selection_frequency(mask_all_ord_freq_selection_run1);
freq_map_run2=compute_selection_frequency(mask_all_ord_freq_selection_run2);

%% test fingerprinting metrics: Idiff, SR, ICC

for rr=1:3 %number of ranges

    TEST_range=squeeze(TEST(:,:,rr));
    RETEST_range=squeeze(RETEST(:,:,rr));
    
    %Identifiability matrix
    Ident_mat=corr(TEST_range,RETEST_range);
    %store for null model testing
    Ident_mat_S_all(:,:,rr) = Ident_mat;
   

    %Success rate: how many times an Iself value is higher than an Iothers
    %value on the same row and column.
    succ_rate = zeros(1,N_subs);
    for ii=1:length(Ident_mat)
        succ_rate(ii) = sum(Ident_mat(ii,ii)>Ident_mat(ii,:))+sum(Ident_mat(ii,ii)>Ident_mat(:,ii));
    end
    sr_indiv_rr_S(:,rr) = (succ_rate*100/(length(Ident_mat)*2-2))'; %individual success rate
    sr_mean = mean(sr_indiv_rr_S(:,rr),2); %average success rate
    
    %Differential identifiability
    Iself = zeros(1,N_subs); %self identifiability
    Iothers = zeros(1,N_subs); % Others-identifiability
    Idiff_indiv = zeros(1,N_subs); %Iself - Iothers

    for s=1:N_subs
        Iself(s) = Ident_mat(s,s);%with mask_diag it was dimension reversed so Stampacchia requires the transposed
        Ident_mat(s,s) = nan;
        Iothers(s) = 0.5.*(nanmean(Ident_mat(s,:)) + nanmean(Ident_mat(:,s))');
        Idiff_indiv(s) = (Iself(s)-Iothers(s));
    end
    Idiff_indiv_rr_S(:,rr)=Idiff_indiv';
    Iself_indiv_rr_S(:,rr)=Iself';
    Iothers_indiv_rr_S(:,rr)=Iothers';
    
   
    %ICC to highlight which roi mainly contributes to the identifiability
    ICC_mat=zeros(N_rois,N_rois);
    ICC = f_ICC_edgewise(TEST_range,RETEST_range);
    ICC_mat(mask_ut) = ICC;
    sum(isnan(ICC))
    any(all(RETEST_range == 0, 2))
    any(all(TEST_range == 0, 2))
    disp('check possible NaN inside ICC'),pause
    ICC_mat = ICC_mat + ICC_mat';
    ICC_mat_reorder=ICC_mat(new_order,new_order);
    ICC_freq(:,:,rr)=ICC_mat_reorder;

    %save for histogram plot
    ICC_nodal_str_hist(:,rr)=sum(ICC_mat); 
    ICC_threshold_reord=prctile(ICC_mat_reorder,icc_mat_thr,'all');

   
    
    %% plot ICC matrix
    
    ICC_mat_S_thr(:,:,rr)=ICC_mat_reorder.*(ICC_mat_reorder>=ICC_threshold_reord);
    ICC_mat_S_full(:,:,rr)=ICC_mat_reorder;

    figure
    imagesc(ICC_mat_S_thr(:,:,rr)); axis square;
    colorbar; 
    w_col=[1 1 1];
    col1=color_range(rr,:);
    col2=[0.5020    0.5020    0.5020];
    [grad1,im]=colorGradient(w_col,col2,200);
    [grad2,im]=colorGradient(w_col,col1,200);
    colormap([grad2])
    title(['ICC Matrix S reconstructed R_' num2str(rr)]);
    hold on
    xticks(new_RSN_bounds+1)
    yticks(new_RSN_bounds+1)
    for jj=1:length(new_RSN_bounds)
        line([new_RSN_bounds(jj)-0.5 new_RSN_bounds(jj)-0.5], [0 79],'Color','k')
        line([0 79],[new_RSN_bounds(jj)-0.5 new_RSN_bounds(jj)-0.5],'Color','k')
    end
    set(gca,'TickLength',[0 0])
    set(gca,'Xtick',[]);set(gca,'Ytick',[]);
   

%     clear S_sub_all_run1 S_sub_all_run2
        
end

%% plot ICC nodal strength distribution

figure
hold on;

for rr = 1:3
    
    %fix number of bins
    x=ICC_nodal_str_hist(:,rr);
    bw=2*iqr(x)/length(x)^(1/3);
    %nbin=round((max(x)-min(x))/bw);
    
    % Plot histogram with transparency
    histogram(x, 'Normalization', 'pdf', ...
        'FaceColor', color_range(rr,:), ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.3, ...
        'BinWidth', bw); % Adjust bin width as needed
    
    % Kernel density estimation
    [f, xi] = ksdensity(x);
    
    % Plot the smooth density curve with fill for shading
    fill(xi, f, color_range(rr,:), ...
        'FaceAlpha', 0.5, ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.2);
end

xlabel('Value');
ylabel('Density');
legend({'Hist 1','KDE 1','Hist 2','KDE 2','Hist 3','KDE 3'});

%% plot Idiff related metrics and SR

%Idiff
figure
DATA_box=[];
for rr=1:3
    DATA_box=[DATA_box;Idiff_indiv_rr_S(:,rr)];
end
group_box=[ones(1,N_subs)';2*ones(1,N_subs)';3*ones(1,N_subs)'];
pos=[1 2 3];
symbol={'o','o','o'};
plot_box_scatter(DATA_box,group_box,pos,color_range,symbol)

title('Differential identifiability')
xlabel('ranges')

%Iself
figure
DATA_box=[];
for rr=1:3
    DATA_box=[DATA_box;Iself_indiv_rr_S(:,rr)];
end
group_box=[ones(1,N_subs)';2*ones(1,N_subs)';3*ones(1,N_subs)'];
pos=[1 2 3];
symbol={'o','o','o'};
plot_box_scatter(DATA_box,group_box,pos,color_range,symbol)

title('Iself')
xlabel('ranges')

%Iothers
figure
DATA_box=[];
for rr=1:3
    DATA_box=[DATA_box;Iothers_indiv_rr_S(:,rr)];
end
group_box=[ones(1,N_subs)';2*ones(1,N_subs)';3*ones(1,N_subs)'];
pos=[1 2 3];
symbol={'o','o','o'};
plot_box_scatter(DATA_box,group_box,pos,color_range,symbol)

title('Iothers')
xlabel('ranges')

%SR
figure
DATA_box=[];
for rr=1:3
    DATA_box=[DATA_box;sr_indiv_rr_S(:,rr)];
end
group_box=[ones(1,N_subs)';2*ones(1,N_subs)';3*ones(1,N_subs)'];
pos=[1 2 3];
symbol={'o','o','o'};
plot_box_scatter(DATA_box,group_box,pos,color_range,symbol)

title('Success rate')
xlabel('ranges')

% disp('check with boxplot function'),pause

%% repeated measures ANOVA for Idiff and SR

%Define repeated-measures factor names
RangeNames = {'Range1', 'Range2', 'Range3'};


% Idiff repeated-measures analysis

% Create table for Idiff
IdiffTable = array2table(Idiff_indiv_rr_S, 'VariableNames', RangeNames);
IdiffTable.Subject = (1:N_subs)';  % subject IDs

% Define repeated measures design
WithinDesign = table(RangeNames', 'VariableNames', {'Range'});
rm_Idiff = fitrm(IdiffTable, 'Range1-Range3 ~ 1', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovatbl_Idiff = ranova(rm_Idiff);

% Post-hoc multiple comparisons for repeated-measures factor 'Range'
c_Idiff = multcompare(rm_Idiff, 'Range');


%Success Rate (sr_indiv_rr_S) analysis

% Create table for Success Rate
srTable = array2table(sr_indiv_rr_S, 'VariableNames', RangeNames);
srTable.Subject = (1:N_subs)';  % subject IDs

% Define repeated measures design
rm_sr = fitrm(srTable, 'Range1-Range3 ~ 1', 'WithinDesign', WithinDesign);

% Run repeated measures ANOVA
ranovatbl_sr = ranova(rm_sr);

% Post-hoc multiple comparisons for repeated-measures factor 'Range'
c_sr = multcompare(rm_sr, 'Range');

%% ICC block matrix visualization (full version)

N_nets_allsub=8;
ICC_block_matrix = zeros(N_nets_allsub, N_nets_allsub, 3);
group_num_allsub=group_num;
group_num_allsub(63:end)=8;

new_order_allsub=[];
for nn=1:N_nets_allsub
    idx_net=find(group_num_allsub==nn);
    new_order_allsub=[new_order_allsub,idx_net];
    
end

RSN_names={'VIS','SOM','DAN','SVAN','LIM','CON','DMN','SUB'};

for rr = 1:3

    ICC_block_matrix(:,:,rr) = computeRSNBlockMatrix(ICC_mat_S_full(:,:,rr), group_num_allsub, new_order_allsub, N_nets_allsub);
    %also compute block frequency DAG matrix for comparison
    freq_block_matrix_run1(:,:,rr) = computeRSNBlockMatrix(freq_map_run1(:,:,rr), group_num_allsub, new_order_allsub, N_nets_allsub);
    freq_block_matrix_run2(:,:,rr) = computeRSNBlockMatrix(freq_map_run2(:,:,rr), group_num_allsub, new_order_allsub, N_nets_allsub);

end

% Plot RSN-level ICC block matrices
for rr = 1:3
    figure('Name', ['ICC Block Matrix Range ', num2str(rr)], 'Color', 'w');
    imagesc(ICC_block_matrix(:,:,rr));
    axis square;
    title(['ICC Block Matrix (Range ', num2str(rr), ')']);
    
    % Custom colormap
    col1 = color_range(rr,:);
    col2 = [0.5020 0.5020 0.5020];
    [grad1,~] = colorGradient([1 1 1], col2, 200);
    [grad2,~] = colorGradient([1 1 1], col1, 200);
    colormap([flipud(grad1); grad2])
    if rr==1
        caxis([-0.1 0.1])
    else
        caxis([-0.4 0.4])
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

%% plot comparison ICC and freq selection

%boxplot
plot_ICC_vs_freq(freq_map_run1, freq_map_run2, ICC_freq, triu_mask);

%block scatter plot
plot_ICC_vs_blockFreq(ICC_block_matrix, freq_block_matrix_run1, freq_block_matrix_run2);


%% Random permutation to test Idiff and SR

%number of null random permutations
null_runs=1000;

for rr=1:3

    Ident_rr=Ident_mat_S_all(:,:,rr);
    
    %create random models with different approaches
    [idiff_median_null,success_median_null, ...
        idiff_mean_null,success_mean_null] = f_ID_permutation_mod_ok(Ident_rr,N_subs,null_runs); 
    
    
    Idiff_perm_all(rr).idiff_median_null=idiff_median_null;
    Idiff_perm_all(rr).idiff_mean_null=idiff_mean_null;

    SR_perm_all(rr).SR_median_null=success_median_null;
    SR_perm_all(rr).SR_mean_null=success_mean_null;
    
   
    %test with median values
    Idiff_median=median(Idiff_indiv_rr_S(:,rr));
    sr_median=median(sr_indiv_rr_S(:,rr));
    thr_idiff_med = prctile(Idiff_perm_all(rr).idiff_median_null, 99);
    results_null_perm(rr).idiff_median_sig = (Idiff_median>thr_idiff_med);
    thr_sr_med = prctile(SR_perm_all(rr).SR_median_null, 99);
    results_null_perm(rr).success_median_sig = (sr_median>thr_sr_med);

    %test with mean values
    Idiff_mean=mean(Idiff_indiv_rr_S(:,rr));
    sr_mean=mean(sr_indiv_rr_S(:,rr));
    thr_idiff_mean = prctile(Idiff_perm_all(rr).idiff_mean_null, 99);
    results_null_perm(rr).idiff_mean_sig = (Idiff_mean>thr_idiff_mean);
    thr_sr_mean = prctile(SR_perm_all(rr).SR_mean_null, 99);
    results_null_perm(rr).success_mean_sig = (sr_mean>thr_sr_mean);


    clear Ident_rr

end


%% comparison with methods of Stampacchia et al., 2024

S_ID_metrics = struct(); % Initialize as a structure
    for rr = 1:3
        [S_ID_metrics(rr).idiff_m, ...
         S_ID_metrics(rr).idiff_se, ...
         S_ID_metrics(rr).idiff_vect, ...
         S_ID_metrics(rr).success, ...
         S_ID_metrics(rr).iself_m, ...
         S_ID_metrics(rr).iself_sd, ...
         S_ID_metrics(rr).Iself_vect, ...
         S_ID_metrics(rr).Iothers_vect, ...
         S_ID_metrics(rr).Iothers_m] = f_ID(Ident_mat_S_all(:,:,rr));
    end

% Null permutation testing
ID_metrics_null = struct(); % Initialize as a structure
    for rr = 1:3
        [ID_metrics_null(rr).idiff_null, ...
         ID_metrics_null(rr).success_null] = f_ID_permutation_orig(Ident_mat_S_all(:,:,rr), null_runs);
        % Test if real IDiff and Success-rate significantly different than chance at 0.01 
        thr_idiff = prctile(ID_metrics_null(rr).idiff_null, 99);
        results_null_perm(rr).idiff_sig = (S_ID_metrics(rr).idiff_m>thr_idiff);
        thr_sr = prctile(ID_metrics_null(rr).success_null, 99);
        results_null_perm(rr).success_sig = (S_ID_metrics(rr).success>thr_sr);
    end


%% surface ICC

clc
clearvars -except ICC_nodal_str_hist
close all

% Convert ICC values to percentile ranks (1-99) for visualization

for rr = 1:3
    
    % Take full nodal ICC
    ICC_values_range = ICC_nodal_str_hist(:, rr);

    % Compute 1st and 99th percentile thresholds
    p1 = prctile(ICC_values_range, 1);
    p99 = prctile(ICC_values_range, 99);

    % Create mask: keep values between 1st and 99th percentiles
    mask = (ICC_values_range >= p1) & (ICC_values_range <= p99);

    % Extract valid values
    valid_vals = ICC_values_range(mask);

    % Compute normalized ranks mapped to 1-99
    percentile_vals = ((tiedrank(valid_vals) - 1) / (numel(valid_vals) - 1)) * 98 + 1;

    % Create array for percentile-transformed values
    ICC_percentile = nan(size(ICC_values_range));
    ICC_percentile(mask) = percentile_vals;

    % Save for later
    ICC_render_strength_rr(:, rr) = ICC_percentile;
end


%% create nifti files for brain renders

%Load atlas file (the atlas order is the same as in EC)
REF_file='/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Schaefer_62Parcels_7Networks_AAL3_FSLMNI152_2mm_2_sym.nii.gz';
atlas_path={REF_file};
  
for rr=1:3
    
    %cortical ROIS
    atlas_Sch = load_untouch_nii(atlas_path{1});
    atlas_Sch = atlas_Sch.img;
    Atlas_templ=zeros(size(atlas_Sch));
    
    rois_cortical=1:62;

    for l= rois_cortical
        
        idx_l=find(atlas_Sch==l);
        Atlas_templ(idx_l)=ICC_render_strength_rr(l,rr);
        
    end
    
    %save nifti file
    save_NII_cortical = create_3D_nii(REF_file,Atlas_templ);
    filename_cort = fullfile('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/FIGURES', ...
        sprintf('ICC_strength_S_R%d_CORTICAL.nii.gz', rr));
    save_untouch_nii(save_NII_cortical,filename_cort);
    
    clear Atlas_templ
    
    %subcortical ROIS
    rois_subcortical=(63:74);
    Atlas_temp_subcortical=zeros(size(atlas_Sch));
    
    for l= rois_subcortical
        
        idx_l_sub=find(atlas_Sch==l);
        Atlas_temp_subcortical(idx_l_sub)=ICC_render_strength_rr(l,rr);
        
    end
    
    %save
    save_NII_subcortical = create_3D_nii(REF_file,Atlas_temp_subcortical);
    filename_sub = fullfile('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/FIGURES', ...
        sprintf('ICC_strength_S_R%d_SUBCORTICAL.nii.gz', rr));
    save_untouch_nii(save_NII_subcortical,filename_sub);
    
    clear Atlas_temp_subcortical

end
