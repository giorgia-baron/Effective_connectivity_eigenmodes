
addpath(genpath('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA'))
addpath(genpath('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/FINAL_VERSION_CODE/functions'))

run1_flag=0;

if run1_flag
    clearvars -except run1_flag
    clc
    close all
    disp('run1')
else
    close all
    clearvars -except run1_flag EIG_EC S_new_decomp E_tot mask_eig_all_sub eigval_all_sub %mean_Sigma_run1  mask_eig_all_run1 eigval_all_run1 Sigma_modes_range_norm_run1 S_modes_range_norm_run1 run1_flag 
    disp('run2')
end

% load data and path to data
load ('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/DATA_EC_HCP.mat')

N_rois=74; %subcorticals included

%% RUN1

%iterate over subjects
for sbj=1:numel(EC_MODES) % -1 se mettiamo soglia su Re
    
    
    clearvars -except run1_flag sbj N_rois EC_MODES eigval_all_sub mask_eig_all_sub S_new_decomp E_tot EIG_EC %NON CANCELLA LE E_mm_real e E_mm_Im dei soggetti!!!!!

    current_sbj=EC_MODES(sbj).reliableID;

    if run1_flag
        EC_sub=EC_MODES(sbj).EC1;
        disp(current_sbj)
    
        %matrix Q
        Q_noiseVar_sub=EC_MODES(sbj).nv1*eye(size(EC_MODES(sbj).EC1));
        
    else
        EC_sub=EC_MODES(sbj).EC2;
        disp(current_sbj)
    
        %matrix Q
        Q_noiseVar_sub=EC_MODES(sbj).nv2*eye(size(EC_MODES(sbj).EC2));

    end

  

    %eigendecomposition of EC_sub
    [V_not_ord,D_not_ord,W_not_ord] = eig(EC_sub,'vector');
    [~,sort_idx_real]=sort(real(D_not_ord),'ascend'); %sort by real part of eigvalues


    %%Sorted eigs
    %eigenvalues
    eigval_sub=D_not_ord(sort_idx_real);
    %right eigenvectors
    selected_V=V_not_ord(:,sort_idx_real);
    

    %evaluate eig couples
    num_couple_eig=couple_eig(eigval_sub,N_rois);
    if sum(num_couple_eig)~=N_rois
        disp('error in computing eigenvalue pairs')
        pause
        break
    end

   

    %compute energy of each mode
    for mm=1:length(num_couple_eig)

        if mm==1
            start_idx=1;
        else
            start_idx=sum(num_couple_eig(1:mm-1))+1;
        end

        range=start_idx:start_idx+num_couple_eig(mm)-1;


        w(mm)=unique(abs(imag(eigval_sub(range)))); 

        gamma(mm)=unique((real(eigval_sub(range)))); 



        E_mm_real(mm)=(gamma(mm)^2)/(-2*gamma(mm));
        E_mm_Im(mm)=(w(mm)^2)/(-2*gamma(mm));



    end


%         E_mm_real_all(sbj).real_energy=E_mm_real;
%         E_mm_Im_all(sbj).im_energy=E_mm_Im;

    %Testing a threshold on Re(lambda) to remove faster modes
    thr=-0.7;
    thr_mask=gamma>thr;
    gamma_thr=gamma(thr_mask);
    E_mm_thr_real=E_mm_real(thr_mask);
    E_mm_thr_Im=E_mm_Im(thr_mask);
    idx_gamma=find(thr_mask,1);
    w_thr=w(thr_mask);

    num_couple_eig_thr=num_couple_eig(thr_mask);

    %%%% Calcolo la AUC considerando la parte REALE degli autovalori

    AUCTOT=trapz(E_mm_real(idx_gamma:end)); 
    ii=0;
    AUC_Re=0;
    while AUC_Re<AUCTOT/2
        AUC_Re=trapz(E_mm_real(idx_gamma:idx_gamma+ii));
        ii=ii+1;
        %AUC_tot(ii)=AUC;
    end
 

    NN_Re=ii-1; %CUTOFF IDX RE
    %NNL=length(num_couple_eig_thr)-NN_Re;


    if NN_Re>0 && NN_Re<length(E_mm_thr_real) 

        cutoff_lambda_sub=gamma_thr(NN_Re);


        % %%% Ordino rispetto alla parte immaginaria positiva crescente
        pos_mask_w=ones(size(w_thr))>0;%w_thr>0; %remove elements with E_Im=0 GB 110625 (verified that it's equal to retain 0 elements in w in terms of cutoff but not of energy) eg isequal([E.cut_im_run2],[E_tot.cut_im_run2]) is ok
        %this is because the first nonzero element of E_mm_thr_Im_sort is
        %very low
        [w_thr_sort,sort_idx_im]=sort(w_thr(pos_mask_w), 'ascend');
        E_mm_thr_posIm=E_mm_thr_Im(pos_mask_w);
        E_mm_thr_Im_sort=E_mm_thr_posIm(sort_idx_im); %ordino i valori di energia immaginaria considerando la parte immaginaria crescente
 
        
        AUCTOT_Im=trapz(E_mm_thr_Im_sort);
        h=0;
        AUC_Im=0;
        while AUC_Im<AUCTOT_Im/2
            AUC_Im=trapz(E_mm_thr_Im_sort(1:1+h));
            h=h+1;
            %AUC_tot_Im(h)=AUC_Im;
        end

        NN_Im=h-1; %CUTOFF IDX IM 

        if NN_Im>0 && NN_Im<length(E_mm_thr_Im_sort)
            cutoff_lambda_sub_Im=w_thr_sort(NN_Im);

            %define ranges of modes & filter subjects
            R_tot=struct();
            R_tot(1).range=gamma_thr<=cutoff_lambda_sub; % corrispondono ai modi con parte reale più BASSA 
            R_tot(2).range=gamma_thr>cutoff_lambda_sub; % corrispondono ai modi con parte reale più ALTA
            R_tot(3).range=w_thr<=cutoff_lambda_sub_Im; % corrispondono ai modi con parte immaginaria più BASSA -> qui sto considerando l'ordine degli autovalori sempre rispetto alla parte Re
            R_tot(4).range=w_thr>cutoff_lambda_sub_Im; % corrispondono ai modi con parte immaginaria più ALTA
            
            %modif GB 250625 change <= with < and viceversa (also
            %below)-->wrong

            % define ranges of modes before removing faster modes (Re<-0.7) ->
            % serve per proiettare correttamente sugli autovalori di
            % interesse andando a mettere zero per tutto gli autovalori al
            % di sotto della soglia sulla parte reale
            
            %now define ranges of original length (=gamma length)
            R_torecon=struct();
            range_re_low=gamma<=cutoff_lambda_sub;
            range_re_low(1:idx_gamma-1)=0;

            range_re_high=gamma>cutoff_lambda_sub;
            range_re_high(1:idx_gamma-1)=0;

            range_im_low=w<=cutoff_lambda_sub_Im;
            range_im_low(1:idx_gamma-1)=0;

            range_im_high=w>cutoff_lambda_sub_Im;
            range_im_high(1:idx_gamma-1)=0;

            range_1=zeros(1,length(range_re_low));
            range_2=zeros(1,length(range_re_low));
            range_3=zeros(1,length(range_re_low));
            range_4=zeros(1,length(range_re_low));

            %Definisco i range combinando i cutoff su Re e Im

            for leng=1:length(range_re_low)

                if range_re_low(leng)==1 && range_im_low(leng)==1
                    range_1(leng)=1;
                end
                if range_re_low(leng)==1 && range_im_high(leng)==1
                    range_2(leng)=1;
                end
                if range_re_high(leng)==1 && range_im_low(leng)==1
                    range_3(leng)=1;
                end
                if range_re_high(leng)==1 && range_im_high(leng)==1
                    range_4(leng)=1;
                end
            end
            
            
            %the following vectors have length=total number of modes as
            %raw unthresholded gamma
            R_torecon(1).range=range_1; % Parte reale minore cutoff e parte immaginaria minore cutoff
            R_torecon(2).range=range_2; % Parte reale minore cutoff e parte immaginaria maggiore cutoff
            R_torecon(3).range=range_3; % Parte reale maggiore cutoff e parte immaginaria minore cutoff
            R_torecon(4).range=range_4; % Parte reale maggiore cutoff e parte immaginaria maggiore cutoff

    
            if sum(R_torecon(1).range)>0 && sum(R_torecon(3).range)>0 && sum(R_torecon(4).range)>0 % mi accerto che tutti i soggetti abbiamo almeno i range 1,3 e 4

                
                %define ranges for eigenvalues
                for r_idx=1:4
                    range_mat=[];
                    range_mask=R_torecon(r_idx).range;

                    %define range for matrix
                    for mm=1:length(num_couple_eig)
                        %range mat=number of rois
                        if range_mask(mm)
                            
                            start_idx = sum(num_couple_eig(1:mm-1)) + 1;
                            end_idx = sum(num_couple_eig(1:mm));
                            range_mat = [range_mat, start_idx:end_idx];

%                             disp('check that range_mat is sequential across ranges, no interruptions')
%                             pause
%                             range_mat=[range_mat,sum(num_couple_eig(1:mm-1))+1:...
%                                 sum(num_couple_eig(1:mm-1))+sum(num_couple_eig(mm))];

                        end

                    end

                    
                    %compute S for each range
                    S_rr=recon_S(selected_V,eigval_sub,Q_noiseVar_sub,range_mat,EC_sub); %rank < number of selected modes
%                     figure,imagesc(imag(S_rr)),colorbar,pause
                    S_real_rr=real(S_rr);
                    S_norm_rr(r_idx)=norm(S_real_rr,'fro');
                    S_modes_range_norm_sub(:,:,r_idx)=S_real_rr./S_norm_rr(r_idx);
                    S_modes_range(:,:,r_idx)=S_real_rr;
               
       
                

                    range_mat_thr=[];
                    %length=number of modes>thr
                    range_mask_thr=R_tot(r_idx).range;

                    %define range to save only eigenvalues of interest
                    %(Re>-0.7) starting from index 1
                    for mm_th=1:length(num_couple_eig_thr)

                        if range_mask_thr(mm_th)

                            start_idx = sum(num_couple_eig_thr(1:mm_th-1)) + 1;
                            end_idx = sum(num_couple_eig_thr(1:mm_th));
                            range_mat_thr = [range_mat_thr, start_idx:end_idx];

%                             range_mat_thr=[range_mat_thr,sum(num_couple_eig_thr(1:mm_th-1))+1:...
%                                 sum(num_couple_eig_thr(1:mm_th-1))+sum(num_couple_eig_thr(mm_th))];

                        end

                    end
                    
%                   imag(eigval_sub(range_mat)), pause
%                   disp('check that range_mat is sequential across ranges,
%                   no interruptions')-->checked, is ok
                    %range_mat_thr=unique(range_mat_thr);
%                     
                  
                 
                    idx_thr_gamma=find(real(eigval_sub)==gamma_thr(1),1);
                    range_mask_tot=zeros(size(eigval_sub(idx_thr_gamma:end)));
                    range_mask_tot(range_mat_thr)=1;
                    if run1_flag
                        mask_eig_all_sub(sbj).mask_eig_run1(r_idx,:)=range_mask_tot; % complex R_tot>thr -0.7
                    else
                        mask_eig_all_sub(sbj).mask_eig_run2(r_idx,:)=range_mask_tot; 
                    end

                    %save dominant eigenmode
                    if r_idx==2 && isempty(range_mat)
                        continue
                    end
                    V_range=selected_V(:,range_mat)';
                    %V_range=V_not_ord(:,range_mat)';
                    D_range=eigval_sub(range_mat);
                    D_range_real=real(D_range);
                    [max_re, idx_max_re]=max(D_range_real)
                    if imag(D_range(idx_max_re(1))<0),disp('error'), pause(),end
                    idx_to_retain=idx_max_re(1);
                    D_range_selected=D_range(idx_to_retain);
                    V_range_filt=V_range(idx_to_retain,:);
                    V_range_dominant(:,r_idx)=V_range_filt;
                    D_range_dominant(r_idx)=D_range_selected;

                end %range
                
                % save eigenvalues
                if run1_flag
                    eigval_all_sub(sbj).eigenall_run1=eigval_sub(idx_thr_gamma:end); %taking eigs with Re>thr
                else
                    eigval_all_sub(sbj).eigenall_run2=eigval_sub(idx_thr_gamma:end); %taking eigs with Re>thr
                end

                if run1_flag
                    S_new_decomp(sbj).IDs=EC_MODES(sbj).reliableID;
                    S_new_decomp(sbj).S_run1=S_modes_range;
                    S_new_decomp(sbj).S_run1_norm=S_modes_range_norm_sub(:,:,:);
                    
                    

                    % Save energy values
                    E_tot(sbj).E_tot_run1=E_mm_thr_real;
                    E_tot(sbj).gamma_tot_run1=gamma_thr;
                    E_tot(sbj).E_tot_run1_Im=E_mm_thr_Im_sort;
                    E_tot(sbj).w_tot_run1=w_thr_sort;
                    E_tot(sbj).cut_re_run1=cutoff_lambda_sub;
                    E_tot(sbj).cut_im_run1=cutoff_lambda_sub_Im;

                    %save dominant modes
                    EIG_EC(sbj).V_range_high_re_run1=V_range_dominant(:,:);
                    EIG_EC(sbj).D_range_all_run1=D_range_dominant(:);

                else
                    S_new_decomp(sbj).IDs=EC_MODES(sbj).reliableID;
                    S_new_decomp(sbj).S_run2=S_modes_range;
                    S_new_decomp(sbj).S_run2_norm=S_modes_range_norm_sub(:,:,:);
                    
                    % Salvo i valori di energia
                    E_tot(sbj).E_tot_run2=E_mm_thr_real;
                    E_tot(sbj).gamma_tot_run2=gamma_thr;
                    E_tot(sbj).E_tot_run2_Im=E_mm_thr_Im_sort;
                    E_tot(sbj).w_tot_run2=w_thr_sort;
                    E_tot(sbj).cut_re_run2=cutoff_lambda_sub;
                    E_tot(sbj).cut_im_run2=cutoff_lambda_sub_Im;

                    %save dominant modes
                    EIG_EC(sbj).V_range_high_re_run2=V_range_dominant(:,:);
                    EIG_EC(sbj).D_range_all_run2=D_range_dominant(:);
                end

            else  %not all ranges of interest included
                if run1_flag

                    S_new_decomp(sbj).IDs=EC_MODES(sbj).reliableID;
                    S_new_decomp(sbj).S_run1=[];
                    S_new_decomp(sbj).S_run1_norm=[];
    
                    eigval_all_sub(sbj).eigenall_run1=[];
    
                    mask_eig_all_sub(sbj).mask_eig_run1=[];

                else

                    S_new_decomp(sbj).IDs=EC_MODES(sbj).reliableID;
                    S_new_decomp(sbj).S_run2=[];
                    S_new_decomp(sbj).S_run2_norm=[];
    
                    eigval_all_sub(sbj).eigenall_run2=[];
    
                    mask_eig_all_sub(sbj).mask_eig_run2=[];
                end
                %continue;
            end
        else
            continue;
        end
    else
        continue;
    end
end


