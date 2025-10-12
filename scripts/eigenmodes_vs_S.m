clc
clear
close all

%load data
load('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/group_scatter.mat');
%including 0 Eim (matlab2022b)
load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/S_Sigma_HCP_110625.mat

addpath(genpath('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/FINAL_VERSION_CODE/functions/'))

N_rois=74; %subcorticals included
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
% new_order=[];
% new_RSN_bounds(1)=1;
% for nn=1:N_nets
%     idx_net=find(group_num==nn);
%     new_order=[new_order,idx_net];
%     new_RSN_bounds(nn+1)=length(new_order)+1;
% end
% color_rois_new_order=color_rois(new_order,:);

%color of each range
color_range=[0    0.4471    0.7412; 0.9294    0.6941    0.1255;0.8510    0.3255    0.0980];


%% Compute correlations
% shift di fase tra regioni i-j e entry della triu Sij
% prodotto dei moduli tra regioni i-j e modulo della entry Sij


for subj=1:numel(EIG_EC)-1

    current_sbj=S_Sigma_new_decomp(subj).IDs;
    disp(current_sbj)

    %% sanity check for range 2 to be removed correctly
    if not(isempty(S_Sigma_new_decomp(subj).S_run1)) && not(isempty(S_Sigma_new_decomp(subj).S_run2))

        S_sub_all_run1=S_Sigma_new_decomp(subj).S_run1;
        V_subj_run1=EIG_EC(subj).V_range_high_re_run1;
        D_subj_run1=EIG_EC(subj).D_range_all_run1;

        if size(S_sub_all_run1,3)==4

            S_sub_all_run1(:,:,2)=[];
            V_subj_run1(:,2)=[];
            D_subj_run1(2)=[];

            if isequal(S_Sigma_new_decomp(subj).S_run1(:,:,1),S_sub_all_run1(:,:,1)) && isequal(S_Sigma_new_decomp(subj).S_run1(:,:,3),S_sub_all_run1(:,:,2)) && isequal(S_Sigma_new_decomp(subj).S_run1(:,:,4),S_sub_all_run1(:,:,3))
                disp('Range 2 removed correctly!')
            end
        end

    else
        continue;
    end

    %% extract S and V phase and mod
    clearvars -except S_sub_all_run1 D_subj_run1 V_subj_run1 subj N_rois current_sbj S_Sigma_new_decomp EIG_EC Corr_single_sub_Rehigh_new_decomp color_range color_rois
    
    if size(V_subj_run1)>0

        for r_idx=1:3

            S_real_subj=S_sub_all_run1(:,:,r_idx);
            S_real_subj=S_real_subj-diag(diag(S_real_subj));
%             if diag(diag(S_real_subj))
%                 disp('S diag is not empty')
%                 pause
%             end

            V_range_subj=squeeze(V_subj_run1(:,r_idx));
            V_range_subj_mod=abs(V_range_subj);
            V_range_subj_angle=rad2deg(angle(V_range_subj));
            if sum(V_range_subj_angle>=360)
                disp('angle>=360')
                pause
            end
            if sum(V_range_subj_angle<0)
                disp('angle<0')
%                 pause
            end

          %initialize upper triangular arrays
          triu_mask=triu(ones(N_rois),1);
%             disp('check triu')
%             pause

          mod_prod=zeros(size(triu_mask));
          S_entry=zeros(size(triu_mask));
          angle_diff=zeros(size(triu_mask));

    
            
            for ii=1:N_rois
                for jj=1:N_rois
                    if triu_mask(ii,jj)==1
                        angle_diff(ii,jj)=V_range_subj_angle(ii)-V_range_subj_angle(jj);
                        mod_prod(ii,jj)=V_range_subj_mod(ii)*V_range_subj_mod(jj);
                        S_entry(ii,jj)=S_real_subj(ii,jj);

                        %convert phase to (-180 180) range
                        if angle_diff(ii,jj)==-180
                            disp('angle=-180')
                            angle_diff(ii,jj)=-angle_diff(ii,jj);
                        end
                       
                        if angle_diff(ii,jj)<-180
                            disp('angle<-180')
                            angle_diff(ii,jj)=angle_diff(ii,jj)+360;
                        end
                         if angle_diff(ii,jj)>180
                            disp('angle>180')
                            angle_diff(ii,jj)=angle_diff(ii,jj)-360;
                         end
                        
                        if or(angle_diff(ii,jj)>=360,angle_diff(ii,jj)<=-360)
                            disp('angle>=360 or angle<=-360')
                            pause
                        end

                    end

                end
            end


            %compute correlations

            angle_diff_triu=angle_diff(triu_mask==1);
            S_entry_triu=S_entry(triu_mask==1);
            [p_corr_shift_entryS,pval_shift_entryS]=corr(angle_diff_triu(:),S_entry_triu(:),'type','Pearson');
%             figure,scatter(angle_diff_triu(:),S_diff_triu(:)),title(num2str(r_idx)),pause,close all

            mod_prod_triu=mod_prod(triu_mask==1);
            S_entry_triu=S_entry(triu_mask==1);
            [p_corr_prodmod_entryS,pval_prodmod_entryS]=corr(mod_prod_triu(:),abs(S_entry_triu(:)),'type','Pearson');
%             figure,scatter(mod_prod_triu(:),abs(S_diff_triu(:))),title(num2str(r_idx)),pause,close all

            

            Corr_single_sub_Rehigh_new_decomp(subj).IDs=current_sbj;
            Corr_single_sub_Rehigh_new_decomp(subj).Sdiff_run1(:,r_idx)=S_entry_triu(:);
            Corr_single_sub_Rehigh_new_decomp(subj).anglediff_run1(:,r_idx)=angle_diff_triu(:);
            Corr_single_sub_Rehigh_new_decomp(subj).modprod_run1(:,r_idx)=mod_prod_triu(:);

            Corr_single_sub_Rehigh_new_decomp(subj).p_corr_shift_entryS_run1(:,r_idx)=p_corr_shift_entryS;
            Corr_single_sub_Rehigh_new_decomp(subj).pval_shift_entryS_run1(:,r_idx)=pval_shift_entryS;

            Corr_single_sub_Rehigh_new_decomp(subj).p_corr_prodmod_entryS_run1(:,r_idx)=p_corr_prodmod_entryS;
            Corr_single_sub_Rehigh_new_decomp(subj).pval_prodmod_entryS_run1(:,r_idx)=pval_prodmod_entryS;

        end
    else
        continue
    end
    
   
end


clearvars -except EIG_EC N_rois color_rois Corr_single_sub_Rehigh_new_decomp color_range color_rois

%% density scatter angle plot


for r_idx=1:3
    X_scat=[];
    Y_scat=[];

    nSubjects = numel(Corr_single_sub_Rehigh_new_decomp);

    %vectorize subjects
    for s=1:nSubjects

        X_scat=[X_scat,Corr_single_sub_Rehigh_new_decomp(s).anglediff_run1(:,r_idx)];
        Y_scat=[Y_scat,Corr_single_sub_Rehigh_new_decomp(s).Sdiff_run1(:,r_idx)];
        
    end

    % subplot(1,3,r_idx)
    
    col1=color_range(r_idx,:);
    grad_range = createGradientMap(col1);  % uses default 200 steps
    cbarLimits=[0 800];
    densityScatterChart_new(X_scat(:),Y_scat(:),50,grad_range,cbarLimits)
    
    %std and mean of corr
    % Preallocate the vectors
    
    vec1 = zeros(nSubjects, 1);
    vec2 = zeros(nSubjects, 1);
    
    for s = 1:nSubjects
        
        vec1(s) = Corr_single_sub_Rehigh_new_decomp(s).p_corr_shift_entryS_run1(r_idx);
        vec2(s) = Corr_single_sub_Rehigh_new_decomp(s).pval_shift_entryS_run1(r_idx);
    end
    
    % Element-wise multiplication
    sig_corr = vec1 .* (vec2<0.05);
    
    % Compute the mean
    mean_corr = mean(sig_corr);

    xlabel('phase shift i-j(v)')
    ylabel('Sij entry')
    xlim([-180 180])
    title(['R_',num2str(r_idx),' mean corr=' num2str(mean_corr)])
    %pause
    
end


%% %% density scatter module plot

close all
% addpath('./functions')
% f=figure
for r_idx=1:3
    X_scat=[];
    Y_scat=[];

    for s=1:nSubjects

        X_scat=[X_scat,Corr_single_sub_Rehigh_new_decomp(s).modprod_run1(:,r_idx)];
        Y_scat=[Y_scat,Corr_single_sub_Rehigh_new_decomp(s).Sdiff_run1(:,r_idx)];
    end
    % subplot(1,3,r_idx)
    col1=color_range(r_idx,:);
    grad_range = createGradientMap(col1);  % uses default 200 steps
    cbarLimits=[0 0.05*10e3];
    densityScatterChart_new(X_scat(:),abs(Y_scat(:)),50,grad_range,cbarLimits)

    %std and mean of corr
    % Preallocate vectors
    
    vec1 = zeros(nSubjects, 1);
    vec2 = zeros(nSubjects, 1);
    
    for s = 1:nSubjects
        vec1(s) = Corr_single_sub_Rehigh_new_decomp(s).p_corr_prodmod_entryS_run1(r_idx);
        vec2(s) = Corr_single_sub_Rehigh_new_decomp(s).pval_prodmod_entryS_run1(r_idx);
    end
    
    % Element-wise multiplication
    sig_corr = vec1 .* (vec2<0.05);
    
    % Compute the mean
    mean_corr = mean(sig_corr);

    
    xlabel('Prod module i*j(v)')
    ylabel('abs(Sij entry)')
    title(['R_',num2str(r_idx),' mean corr=' num2str(mean_corr)])
%     pause
    
end


%% plot group-level dominant eigenvector in polar coordinates

for rr=1:3

    n_points = length(squeeze(EIG_EC(1).V_range_high_re_run1(:,1))); % example guess
    n_subjs = numel(EIG_EC)-1;
    V_subj_shift = zeros(n_points, n_subjs);


    %remove mean from each subj and compute group-level V
    for ss=1:numel(EIG_EC)-1
        
        V_subj_run1=EIG_EC(ss).V_range_high_re_run1;
        V_subj_run1(:,2)=[];
        V_subj=squeeze(V_subj_run1(:,rr));
        

        theta_first=angle(mean(V_subj)); %to apply rotation
        V_subj_shift(:,ss)=V_subj*exp(-1i * (theta_first)); %negative sign added to anti rotate

    end


    V_range_mean=mean(V_subj_shift,2);

 
    figure(rr)
    c=compass(V_range_mean)
%      axis([0 0.3 0 0.3]);
    for aa=1:N_rois
        c(aa).Color=color_rois(aa,:);
        c(aa).LineWidth=1.5;
    end
    title(['ROI eigs, range=' num2str(rr)])
    hold off
%     pause()
%     close all
end