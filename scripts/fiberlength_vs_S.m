clear all
clc
close all

%load data
load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/group_scatter.mat');
load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/S_Sigma_HCP_110625.mat
addpath(genpath('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/FINAL_VERSION_CODE/functions/'))

N_rois=74; %subcorticals included
N_subs=143;
N_nets=13;


%get ROI colors
group_col = [0.494117647058824   0.184313725490196   0.556862745098039;0.301960784313725   0.745098039215686   0.933333333333333;...
    0.466666666666667   0.674509803921569   0.188235294117647;1     0     1;0.909803921568627   0.909803921568627   0.619607843137255;...
    0.929411764705882   0.694117647058824   0.125490196078431;0.850980392156863   0.325490196078431   0.098039215686275;1 0 0.5; 0 0 0.5;...
    1 0 0;0.75 0.75 0;0 1 0.5;1.0000    1.0000    0.0667];

color_rois=[];
for rr=1:N_rois
    color_rois(rr,:)=group_col(group_num(rr),:);
end

%reorder labels and colors
new_order=[];
new_RSN_bounds(1)=1;
for nn=1:N_nets
    idx_net=find(group_num==nn);
    new_order=[new_order,idx_net];
    new_RSN_bounds(nn+1)=length(new_order)+1;
end
color_rois_new_order=color_rois(new_order,:);

%color of each range
color_range=[0    0.4471    0.7412; 0.9294    0.6941    0.1255;0.8510    0.3255    0.0980];

%% binarize DAGs and keep reliable links

%load S derived DAGs (RUN1)
dag_r1=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_RUN1_range1.mat').DAG_all_range1;
dag_r2=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_RUN1_range2.mat').DAG_all_range2;
dag_r3=load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/Null_model/RUN1/DAG_data_RUN1_range3.mat').DAG_all_range3;



for range=1:3

    for ss=1:N_subs
        
    
        %binary DAG R1
        DAG_bin=not(eval(['dag_r',num2str(range)]).DAG(:,:,:,ss)==0); 
        DAG_bin_sum=sum(DAG_bin,3);
       
        %keep only links present in all iterations
        tmp=DAG_bin_sum==100;
        %transpose matrix
        tmp=tmp';
        mask_DAG(:,:,ss,range)=tmp;

      
    
    end

end

%% load structural data

for ss=1:N_subs
    
    %fiber length
    fl_path='/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/HCP_young_4_DCM/fiber_length';
    ID_curr=S_Sigma_new_decomp(ss).IDs;
    if isfile(fullfile(fl_path,[ID_curr,'.Schaefer_AAL3_4_DCM_EdgeLength.csv']))

        fl_curr=readmatrix(fullfile(fl_path,[ID_curr,'.Schaefer_AAL3_4_DCM_EdgeLength.csv']));

        %keep ROI order consistency with S
        fl_curr_ord=fl_curr([1:66,68,67,69:74],[1:66,68,67,69:74]);
        SC_FL(:,:,ss)=fl_curr_ord; 
        SC_FL_r1(:,:,ss)=fl_curr_ord.*mask_DAG(:,:,ss,1);
        SC_FL_r2(:,:,ss)=fl_curr_ord.*mask_DAG(:,:,ss,2);
        SC_FL_r3(:,:,ss)=fl_curr_ord.*mask_DAG(:,:,ss,3);
        
    else
        
        SC_FL(:,:,ss)=nan*ones(N_rois);
        SC_FL_r1(:,:,ss)=nan*ones(N_rois);
        SC_FL_r2(:,:,ss)=nan*ones(N_rois);
        SC_FL_r3(:,:,ss)=nan*ones(N_rois);
        

    end
    
    %do the same with number of streamlines
    ns_path='/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/HCP_young_4_DCM/SC';
    if isfile(fullfile(ns_path,[ID_curr,'.Schaefer_AAL3_4_DCM_SC.csv']))

        ns_curr=readmatrix(fullfile(ns_path,[ID_curr,'.Schaefer_AAL3_4_DCM_SC.csv']));
        ns_curr_ord=ns_curr([1:66,68,67,69:74],[1:66,68,67,69:74]);
        SC_NS(:,:,ss)=ns_curr_ord; 
        
    else
        ss
        SC_NS(:,:,ss)=nan*ones(N_rois);
       
    end
    
end

%% mask to keep only reliable structural links and compute group thresholds

% keep links present in at least 90% of subjects
thr=(N_subs-1)*0.9; %one subject has no SC
SC_NS_stable=sum((SC_NS>0),3,'omitnan');
SC_NS_stable_mask=SC_NS_stable>=thr;

%remove links of cerebellum with cortex
SC_NS_stable_mask(67,1:62)=0;
SC_NS_stable_mask(1:62,67)=0;

SC_NS_stable_mask(74,1:62)=0;
SC_NS_stable_mask(1:62,74)=0;

figure,imagesc(SC_NS_stable_mask(new_order,new_order))
hold on
xticks(new_RSN_bounds+1)
yticks(new_RSN_bounds+1)
for jj=1:length(new_RSN_bounds)
    line([new_RSN_bounds(jj)-0.5 new_RSN_bounds(jj)-0.5], [0 75],'Color','k')
    line([0 75],[new_RSN_bounds(jj)-0.5 new_RSN_bounds(jj)-0.5],'Color','k')
end

%get threshold to define FL ranges
SC_FL_mean=mean(SC_FL.*SC_NS_stable_mask,3,'omitnan');
q1_thr_fl_mean=prctile(SC_FL_mean(SC_NS_stable_mask>0),25);
q4_thr_fl_mean=prctile(SC_FL_mean(SC_NS_stable_mask>0),75);

% figure,histogram(SC_FL_mean(SC_NS_stable_mask>0),100),hold on,xline(q1_thr_fl_mean),hold on,xline(q4_thr_fl_mean),pause(),title('FL dist')

%% FL distribution
% Original histogram

figure;
hold on;

% Extract data
data = SC_FL_mean(SC_NS_stable_mask > 0);

% Plot histogram with transparency
histogram(data, 'Normalization', 'pdf', ...
    'FaceColor', [0.2 0.6 0.8], ...   % Set your desired color
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.3, ...
    'BinWidth', 5);               % Adjust bin width as needed

% Kernel density estimation
[x_density, y_density] = ksdensity(data);

% Plot the smooth density curve with fill for shading
hold on
fill(y_density, x_density, [0.2 0.6 0.8], ...
    'FaceAlpha', 0.5, ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.2);

% Add threshold lines
hold on
xline(q1_thr_fl_mean, 'k--', 'LineWidth', 1.5);
xline(q4_thr_fl_mean, 'k--', 'LineWidth', 1.5);

% Labels and title
xlabel('Value');
ylabel('Density');
title('FL dist');

% Legend
legend({'Histogram', 'Smoothed Density', 'Q1 Threshold', 'Q4 Threshold'});

%% compute subject-level intervals

for ss=1:N_subs

    
    if ss==33 %no SC data
        continue
    else

        %get the FL class based on group-based prct
        fl_curr=SC_FL(:,:,ss).*SC_NS_stable_mask;
        fl_short=logical((fl_curr<q1_thr_fl_mean).*SC_NS_stable_mask);
        fl_med=and(fl_curr>=q1_thr_fl_mean,fl_curr<=q4_thr_fl_mean);
        fl_long=fl_curr>q4_thr_fl_mean;
       
        
         %FL (need to mask to compute num_tot correctly)
         SC_FL_r1(:,:,ss)=SC_FL_r1(:,:,ss).*SC_NS_stable_mask;
         SC_FL_r2(:,:,ss)=SC_FL_r2(:,:,ss).*SC_NS_stable_mask;
         SC_FL_r3(:,:,ss)=SC_FL_r3(:,:,ss).*SC_NS_stable_mask;
    
    
         num_tot_fl_r1=sum(SC_FL_r1(:,:,ss)>0,'all');
         num_tot_fl_r2=sum(SC_FL_r2(:,:,ss)>0,'all');
         num_tot_fl_r3=sum(SC_FL_r3(:,:,ss)>0,'all');
    
    
         %count percentage of included links in DAG for each FL class
    
         FL_num_r1_short(ss)=sum(SC_FL_r1(:,:,ss).*fl_short>0,'all')/num_tot_fl_r1;
         FL_num_r2_short(ss)=sum(SC_FL_r2(:,:,ss).*fl_short>0,'all')/num_tot_fl_r2;
         FL_num_r3_short(ss)=sum(SC_FL_r3(:,:,ss).*fl_short>0,'all')/num_tot_fl_r3;
    
         FL_num_r1_med(ss)=sum(SC_FL_r1(:,:,ss).*fl_med>0,'all')/num_tot_fl_r1;
         FL_num_r2_med(ss)=sum(SC_FL_r2(:,:,ss).*fl_med>0,'all')/num_tot_fl_r2;
         FL_num_r3_med(ss)=sum(SC_FL_r3(:,:,ss).*fl_med>0,'all')/num_tot_fl_r3;
         
         FL_num_r1_long(ss)=sum(SC_FL_r1(:,:,ss).*fl_long>0,'all')/num_tot_fl_r1;
         FL_num_r2_long(ss)=sum(SC_FL_r2(:,:,ss).*fl_long>0,'all')/num_tot_fl_r2;
         FL_num_r3_long(ss)=sum(SC_FL_r3(:,:,ss).*fl_long>0,'all')/num_tot_fl_r3;
    
         %sanity check
         if not(FL_num_r1_short+FL_num_r1_med+FL_num_r1_long==1)
             pause
         end
    
         if not(FL_num_r2_short+FL_num_r2_med+FL_num_r2_long==1)
             pause
         end
    
         if not(FL_num_r3_short+FL_num_r3_med+FL_num_r3_long==1)
             pause
         end

    end


end

%% Stacked barplot

% data = [FL_num_r1_short;FL_num_r1_med;FL_num_r1_long]';  % Transpose to get 143x3 matrix
% %remove subj 33
% data(33,:)=[];
% 
% % Create the stacked bar plot
% figure;
% bar(data, 'stacked');
% 
% % Customize appearance
% xlabel('Subject');
% ylabel('Value');
% title('Stacked Bar Plot by Subject-FL');
% legend({'short', 'med', 'long'}, 'Location', 'northeastoutside');
% 
% % Optional: improve readability if too dense
% xlim([0 144]);
% set(gca, 'XTick', 1:10:143);  % Adjust x-axis ticks for readability
% box off;


%% raincloud plots across ranges
close all

%select FL class
sc_class=3

switch sc_class
    case 1 %short
        
            rain_data=[FL_num_r1_short;FL_num_r2_short;FL_num_r3_short];   
        
    case 2 %med
        
            rain_data=[FL_num_r1_med;FL_num_r2_med;FL_num_r3_med];   
          
    case 3 %long
         
            rain_data=[FL_num_r1_long;FL_num_r2_long;FL_num_r3_long];   
       
end

%remove subj 33 because of nan values
rain_data(:,33)=[];
raincloud_plot((rain_data'))
% raincloud_plot(log(rain_data'))


%% statistical testing


%statistical testing (Friedman+post hoc test)
% figure,kruskal_dunn_test(rain_data')
pval_thr=0.01;
friedman_posthoc_test(rain_data',pval_thr)

