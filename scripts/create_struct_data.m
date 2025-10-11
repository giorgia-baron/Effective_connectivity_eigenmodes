clear all
clc
close all

%% load data RUN1

files_run1=dir('/nfsd/biopetmri4/Users/Giorgia/HCP_project/RESULTS/RESULTS_YA_new_1000_RUN1_noisevar');
files_run1=files_run1(3:end);

files_run2=dir('/nfsd/biopetmri4/Users/Giorgia/HCP_project/RESULTS/RESULTS_YA_new_1000_RUN2');
files_run2=files_run2(3:end);

%check if noise variance is right

subj=0;
for subj0=1:length(files_run1)
    subj0
    load(fullfile(files_run1(subj0).folder,files_run1(subj0).name),'A_hat','NoiseVar','y','tmask','sys_est')

    %check noise variance
    noiseVar=diag(sys_est.R)./var(y)';
    if not(noiseVar==0.1*ones(74,1))
        subj
        disp('error')
        pause
    end
   
end

subj=0;
for subj0=1:length(files_run2)
    subj0
    load(fullfile(files_run2(subj0).folder,files_run2(subj0).name),'A_hat','NoiseVar','y','tmask','sys_est')

    %check noise variance
    noiseVar=diag(sys_est.R)./var(y)';
    if not(noiseVar==0.1*ones(74,1))
        subj
        disp('error')
        pause
    end
   
end


%% save reliable data for the two HCP runs

clear all
clc

folder_RUN1='/nfsd/biopetmri4/Users/Giorgia/HCP_project/RESULTS/RESULTS_YA_new_1000_RUN1_noisevar';
folder_RUN2='/nfsd/biopetmri4/Users/Giorgia/HCP_project/RESULTS/RESULTS_YA_new_1000_RUN2';

RUN1=dir(folder_RUN1);
RUN1=RUN1(3:end);

RUN2=dir(folder_RUN2);
RUN2=RUN2(3:end);

subj_ID=1;
discard_subj=[];

for subj=1:length(RUN1) %iterate over smallest folder
    
    RUN1(subj).name

    %check that subj is present in RUN2 folder
    if not(isfile(fullfile(folder_RUN2,RUN1(subj).name)))
        discard_subj=[discard_subj;RUN1(subj).name];
        disp('only run 1 is present')
        continue
    
    else

        %both runs are present
       
        %RUN1
        load(fullfile(folder_RUN1,RUN1(subj).name),'A','tmask','num_iter','NoiseVar')
        EC1=A;
        noiseVar1=NoiseVar;
        if sum(tmask)<432 || num_iter<10 %at least (more or less) 50% of volumes 
            disp('unstable results')
            continue
        end
        
        %RUN2
        load(fullfile(folder_RUN2,RUN1(subj).name),'A','tmask','num_iter','NoiseVar')
        EC2=A;
        noiseVar2=NoiseVar;
        if sum(tmask)<432 || num_iter<10 %at least (more or less) 50% of volumes
            disp('unstable results')
            continue 
        end
        
        if not(isempty(find(EC1))) && not(isempty(find(EC2))) %check they are not equal to zero
            subj_ID
        
            ID=RUN1(subj).name;
            ID=erase(ID,'EMmulti_HCP_');
            ID=erase(ID,'.mat');
            EC_MODES(subj_ID).reliableID=ID;
            EC_MODES(subj_ID).EC1=EC1;
            EC_MODES(subj_ID).nv1=noiseVar1;
            EC_MODES(subj_ID).EC2=EC2;
            EC_MODES(subj_ID).nv2=noiseVar2;

            subj_ID=subj_ID+1;

        end
        
        clear ID EC1 EC2 noiseVar1 noiseVar2 num_iter tmask A

    end

end
