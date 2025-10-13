clear
clc
close all

%% setting params
%only cortical ROIs, otherwise N_nodes=1:74 to include subcortical nodes
N_nodes=1:74;
n=max(N_nodes);

N_subs=143;

addpath('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES');

load ('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/S_Sigma_HCP.mat')
load('/nfsd/biopetmri4/Users/ClaudiaTarricone/EC_modes_hierarchies/Data/node_labels_Schaefer.mat')


WR_tot_range3 = [];
WR_rand_tot_range3 = [];

%% find cycles and sort DAG
range=3;
for ss=1:N_subs
    ss

    S_sub_all=S_Sigma_new_decomp(ss).S_run1;
    S_sub_all(:,:,2)=[];

    if (S_Sigma_new_decomp(ss).S_run1(:,:,1)==S_sub_all(:,:,1)) & (S_Sigma_new_decomp(ss).S_run1(:,:,3)==S_sub_all(:,:,2)) & (S_Sigma_new_decomp(ss).S_run1(:,:,4)==S_sub_all(:,:,3))
        disp('Range 2 removed correctly!')
    end

    S_range=S_sub_all(:,:,range);
    S_sub=S_range'; %matlab convention: outgoing over rows and incoming over columns, need to transpose -> nel loro codice non facevano il trasposto
    S_sub(S_sub<0)=0; 
    S_sub(logical(eye(size(S_sub)))) = 0;
    S_sub_rand = S_sub;
    for i=1:n
        for j=1:i
            if S_sub(i,j)~=0 || S_sub(j,i)~=0
                if rand > 0.5
                    S_sub_rand(i,j) = S_sub(j,i); S_sub_rand(j,i) = S_sub(i,j);
                end
            end
        end
    end

    MATRIX_S_range3(:,:,ss)=S_sub;
    MATRIX_S_rand_range3(:,:,ss)=S_sub_rand;
    all_edges_S=sum(length(nonzeros(S_sub)));
    all_weights_S=sum(nonzeros(S_sub));
    S_sub=digraph(abs(S_sub),node_labels); %create directed graph

    all_edges_S_rand=sum(length(nonzeros(S_sub_rand)));
    all_weights_S_rand=sum(nonzeros(S_sub_rand));
    S_sub_rand=digraph(abs(S_sub_rand),node_labels); %create directed graph
    

    for iter=1:100
        iter
        %set seed
        rng(iter,'Threefry')

        all_edges=all_edges_S;
        all_weights=all_weights_S;
        A_sub=S_sub;
        edges_removed=0;
        counter=1;
        
        while(hascycles(A_sub))
                
               [cycles,edgecycles] = allcycles(A_sub,'MaxNumCycles',100);
               idx_cycle_rand=randperm(length(edgecycles),1); %take random cycle and remove the lowest link
               w_idx=cell2mat(edgecycles(idx_cycle_rand,:));
               min_idx=1;
               w_min=A_sub.Edges.Weight(w_idx(min_idx));
               for ll=2:length(w_idx)
                      w_test=A_sub.Edges.Weight(w_idx(ll));
                       if abs(w_test)<abs(w_min)
                           w_min=w_test;
                           min_idx=ll;
                       end
               end
              
               cycles_rmd_subj(counter)=w_min;
               nodes_cycles_rmd_subj(counter,:)=A_sub.Edges(w_idx(min_idx),:).EndNodes;
               edges_removed=edges_removed+1;
               counter=counter+1;
               A_sub=rmedge(A_sub,w_idx(min_idx));
               clear edgecycles
        end

        if hascycles(S_sub)
            cycles_rmd_range3(iter,ss).weights=cycles_rmd_subj;
            cycles_rmd_range3(iter,ss).nodes=nodes_cycles_rmd_subj;
            WR_perc_range3(iter,ss)=sum(cycles_rmd_subj)/all_weights*100;
            WR_range3(iter,ss)=sum(cycles_rmd_subj);
            ER_perc_range3(iter,ss)=edges_removed/all_edges*100;
            ER_range3(iter,ss)=edges_removed;
        else
            cycles_rmd_range3(iter,ss).weights=[];
            cycles_rmd_range3(iter,ss).nodes=[];
            WR_perc_range3(iter,ss)=NaN;
            WR_range3(iter,ss)=NaN;
            ER_perc_range3(iter,ss)=NaN;
            ER_range3(iter,ss)=NaN;
        end
    
%         cycles_rmd_range3(iter,ss).weights=cycles_rmd_subj;
%         cycles_rmd_range3(iter,ss).nodes=nodes_cycles_rmd_subj;
         
        %save DAG matrix
        A_sub_mat_range3=full(adjacency(A_sub,"weighted"));
        DAG_all_range3.DAG(:,:,iter,ss)=A_sub_mat_range3; %transposed wrt EC
       
%         %save removed edges
%         WR_perc_range3(iter,ss)=sum(cycles_rmd_subj)/all_weights*100;
%         WR_range3(iter,ss)=sum(cycles_rmd_subj);
%         ER_perc_range3(iter,ss)=edges_removed/all_edges*100;
%         ER_range3(iter,ss)=edges_removed;
    
        clear cycles_rmd_subj nodes_cycles_rmd_subj A_sub_mat_rand_range3 w_idx w_min A_sub all_edges all_weights edges_removed counter


        all_edges_rand=all_edges_S_rand;
        all_weights_rand=all_weights_S_rand;
        A_sub_rand=S_sub_rand;
        edges_removed=0;
        counter=1;
        
        while(hascycles(A_sub_rand))
                
               [cycles,edgecycles] = allcycles(A_sub_rand,'MaxNumCycles',100);
               idx_cycle_rand=randperm(length(edgecycles),1); %take random cycle and remove the lowest link
               w_idx=cell2mat(edgecycles(idx_cycle_rand,:));
               min_idx=1;
               w_min=A_sub_rand.Edges.Weight(w_idx(min_idx));
               for ll=2:length(w_idx)
                      w_test=A_sub_rand.Edges.Weight(w_idx(ll));
                       if abs(w_test)<abs(w_min)
                           w_min=w_test;
                           min_idx=ll;
                       end
               end
              
               cycles_rmd_subj(counter)=w_min;
               nodes_cycles_rmd_subj(counter,:)=A_sub_rand.Edges(w_idx(min_idx),:).EndNodes;
               edges_removed=edges_removed+1;
               counter=counter+1;
               A_sub_rand=rmedge(A_sub_rand,w_idx(min_idx));
               clear edgecycles
        end
        if hascycles(S_sub_rand)
                 cycles_rmd_rand_range3(iter,ss).weights=cycles_rmd_subj;
                 cycles_rmd_rand_range3(iter,ss).nodes=nodes_cycles_rmd_subj;
                 WR_perc_rand_range3(iter,ss)=sum(cycles_rmd_subj)/all_weights_rand*100;
                 WR_rand_range3(iter,ss)=sum(cycles_rmd_subj);
                 ER_perc_rand_range3(iter,ss)=edges_removed/all_edges_rand*100;
                 ER_rand_range3(iter,ss)=edges_removed;
        else
                 cycles_rmd_rand_range3(iter,ss).weights=[];
                 cycles_rmd_rand_range3(iter,ss).nodes=[];
                 WR_perc_rand_range3(iter,ss)=NaN;
                 WR_rand_range3(iter,ss)=NaN;
                 ER_perc_rand_range3(iter,ss)=NaN;
                 ER_rand_range3(iter,ss)=NaN;
        end
%         cycles_rmd_rand_range3(iter,ss).weights=cycles_rmd_subj;
%         cycles_rmd_rand_range3(iter,ss).nodes=nodes_cycles_rmd_subj;
%         
        %save DAG matrix
        A_sub_mat_rand_range3=full(adjacency(A_sub_rand,"weighted"));
        DAG_all_rand_range3.DAG(:,:,iter,ss)=A_sub_mat_rand_range3; %transposed wrt EC
       
%         %save removed edges
%         WR_perc_rand_range3(iter,ss)=sum(cycles_rmd_subj)/all_weights_rand*100; 
%         WR_rand_range3(iter,ss)=sum(cycles_rmd_subj);
%         ER_perc_rand_range3(iter,ss)=edges_removed/all_edges_rand*100;
%         ER_rand_range3(iter,ss)=edges_removed;
    
        clear cycles_rmd_subj nodes_cycles_rmd_subj A_sub_mat_rand_range3 w_idx w_min A_sub_rand all_edges_rand all_weights_rand edges_removed counter

    end

    WR_tot_range3 = [WR_tot_range3 min(WR_range3(:,ss))];
    WR_rand_tot_range3 = [WR_rand_tot_range3 min(WR_rand_range3(:,ss))]; 

end

out_path=('/nfsd/biopetmri4/Users/ClaudiaTarricone/EC_modes_hierarchies/Results/Null_model/RUN1');
save(fullfile(out_path,'DAG_data_RUN1_range3.mat'),'DAG_all_range3','MATRIX_S_range3');
save(fullfile(out_path,'Results_test_null_DAG_RUN1_range3.mat'),'WR_tot_range3','WR_rand_tot_range3','WR_rand_range3','WR_range3','ER_perc_rand_range3','ER_rand_range3','ER_range3','ER_perc_range3','cycles_rmd_range3','cycles_rmd_rand_range3')

%%
f=figure;
boxplot([WR_tot_range3',WR_rand_tot_range3'],'Labels',{'original','randomized'})
title(['Weights removed-run1(143 subjects)-range ' num2str(range)])
set(f, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

path_save=('/nfsd/biopetmri4/Users/ClaudiaTarricone/EC_modes_hierarchies/Results/Null_model/RUN1');

saveas(f,fullfile(path_save,['Comparison_removed_edges_null_model_RUN1_range3.png'])) 
%save('...','DAG','MATRIX_EC_ANTISYM','MATRIX_S',...
%    'ER','ER_perc','cycles_rmd')
