%% save DCM fit

% clear all
% clc
% close all
% %% 
% load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/DATA_EC_HCP.mat
% 
% subj=0;
% folder = '/nfsd/biopetmri4/Users/Giorgia/HCP_project/RESULTS/RESULTS_YA_new_1000_RUN2/'
% 
% %% RUN 1
% 
% for sub=1:length(EC_MODES)
% 
%    sub
%    
%    subj_name=EC_MODES(sub).reliableID;
%    data_path=fullfile(folder,['EMmulti_HCP_' subj_name '.mat']);
%    if not(isfile(data_path))
%        disp('error')
%        pause
%    end
%   
%     
%    subj=subj+1;
%    
%    %load corr eFC-sFC
%    corr_FC(subj)=load(data_path).FC_p_corr;
% 
%    %load y
%    y=load(data_path).y;
%    y_pred=load(data_path).y_pred; %Kalman one-step prediction
%    tmask=load(data_path).tmask;
%     
%    
%    Y(:,:,subj)=y;
%    Y_PRED(:,:,subj)=y_pred;
%    NAME{subj}=subj_name;
%    RMSE(:,subj)=rmse(y_pred,y);
%    TMASK(:,subj)=tmask;
%     
%   
% end

%% create boxchart for RMSE

clc
close all
clear all

%load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/DCM_fit_RUN1.mat
load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/DCM_fit_RUN2.mat

load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/node_labels_Schaefer.mat

% new_idx = [36:37, 72:74, 30:35, 65:71, 25:29, 61:64, 21:24, 55:60, 18:20, 53:54, 13:17, 45:52, 1:12, 38:44, 75, 81, 76,82, 77, 83, 78, 84, 79, 86, 80, 85]; %vis moto dors sal limbico cont dmn sub (Thalamus, Caudate, Putamen, Pallidum, Cerebellum, Hippocampus) 
% node_labels=node_labels(new_idx);
% width = [5, 13, 9, 10, 5, 13, 19, 2, 2, 2, 2, 2, 2];
% M0 =[ones(5,1); 2*ones(13,1); 3*ones(9,1); 4*ones(10,1); 5*ones(5,1); 6*ones(13,1); 7*ones(19,1); 8*ones(12,1)];
% lines = cumsum(width)+0.5;
% lines=lines(1:end-1);

%compute NRMSE
for ss=1:length(NAME)

    y=Y(:,:,ss);
    y_pred=Y_PRED(:,:,ss);
    y_min=prctile(y,25);
    y_max=prctile(y,75);
    NRMSE(:,ss)=(RMSE(:,ss)./(y_max-y_min)')*100;

end

% ordina le roi
% NRMSE = NRMSE(new_idx,:);


figure(1)
boxplot(NRMSE','Colors','k','Symbol','')
ylabel('NRMSE')
title('Dataset 2')
% xline(lines)
% xlabel('ROIs')

hold on

ROI_node=1:74;

for rr=1:max(ROI_node)
    
    rr
    for ss=1:length(NAME)
        
        a=-0.2;
        b=0.2;
        rn=(b-a).*rand + a;
        
        scatter(ROI_node(rr)+rn,NRMSE(rr,ss),12,'r','filled')
        hold on
    
    end

end

median_FC = median(corr_FC);
median_NRMSE = median(NRMSE(:));
% median_NRMSE_sbj = median(NRMSE);
% median_NRMSE_tot = median(median_NRMSE_sbj);

%% create boxchart for FCcorr

figure(2)
boxplot(corr_FC,'Symbol','','Labels',{'Dataset 2'})
ylabel('corr eFC-sFC')
hold on

for ss=1:length(corr_FC)
    
    a=-0.1;
    b=0.1;
    rn=(b-a).*rand + a;
   
    scatter(1+rn,corr_FC(ss),15,'r','filled')
    


end

   


