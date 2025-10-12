clear
clc
close all

%% setting params
%only cortical ROIs, otherwise N_nodes=1:74 to include subcortical nodes
N_nodes=1:74;
n=max(N_nodes);

% node_labels=node_labels(N_nodes);
N_subs=144;
% EC_all=EC_all(N_nodes,N_nodes,:);
% group_num=group_num(N_nodes);
% group=group(N_nodes);

%files=dir('./DATA/example_data_4_DAG'); %path to data
%files=files(3:end);

load ('/nfsd/biopetmri4/Users/ClaudiaTarricone/EC_modes_hierarchies/Data/node_labels_Schaefer.mat');
load ('/nfsd/biopetmri4/Users/ClaudiaTarricone/EC_modes_hierarchies/Data/all_S.mat');


WR_tot = [];
WR_rand_tot = [];

%% find cycles and sort DAG

for ss=1:N_subs
    ss
    
%     EC_sub=load(fullfile(files(ss).folder,files(ss).name)).A; 
%     EC_sub=EC_sub(N_nodes,N_nodes);
% 
%     % estimate matrix S
%     NoiseVar2=load(fullfile(files(ss).folder,files(ss).name)).NoiseVar;
%     Q_sub=NoiseVar2.*eye(n);
%     Sigma_sub=lyap(EC_sub,Q_sub); %zero-lag covariance matrix
%     S_sub=(EC_sub*Sigma_sub-Sigma_sub*EC_sub'); %time-lag covariance matrix S
%     EC_sub_check=0.5*(-Q_sub+S_sub)*inv(Sigma_sub); %test if EC_sub_check=EC_sub
%     S_sub=S_sub'; %matlab convention: outgoing over rows and incoming over columns, need to transpose
%     S_sub(S_sub<0)=0; % set negative weights to zero
%     MATRIX_S(:,:,ss)=S_sub;

    S_sub = S_all_subj(ss).S_RUN1;
    S_sub(S_sub<0)=0;
    
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

    MATRIX_S(:,:,ss)=S_sub;
    MATRIX_S_rand(:,:,ss)=S_sub_rand;
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
    
        cycles_rmd(iter,ss).weights=cycles_rmd_subj;
        cycles_rmd(iter,ss).nodes=nodes_cycles_rmd_subj;
        
        %save DAG matrix
        A_sub_mat=full(adjacency(A_sub,"weighted"));
        DAG(:,:,iter,ss)=A_sub_mat; %transposed wrt EC
       
        %save removed edges
        WR_perc(iter,ss)=sum(cycles_rmd_subj)/all_weights*100;
        WR(iter,ss)=sum(cycles_rmd_subj);
        ER_perc(iter,ss)=edges_removed/all_edges*100;
        ER(iter,ss)=edges_removed;
    
        clear cycles_rmd_subj nodes_cycles_rmd_subj A_sub_mat w_idx w_min A_sub all_edges all_weights edges_removed counter


        all_edges_rand=all_edges_S_rand;
        all_weights_rand=all_weights_S_rand;
        A_sub=S_sub_rand;
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
    
        cycles_rmd_rand(iter,ss).weights=cycles_rmd_subj;
        cycles_rmd_rand(iter,ss).nodes=nodes_cycles_rmd_subj;
        
        %save DAG matrix
        A_sub_mat=full(adjacency(A_sub,"weighted"));
        DAG_rand(:,:,iter,ss)=A_sub_mat; %transposed wrt EC
       
        %save removed edges
        WR_perc_rand(iter,ss)=sum(cycles_rmd_subj)/all_weights_rand*100; 
        WR_rand(iter,ss)=sum(cycles_rmd_subj);
        ER_perc_rand(iter,ss)=edges_removed/all_edges_rand*100;
        ER_rand(iter,ss)=edges_removed;
    
        clear cycles_rmd_subj nodes_cycles_rmd_subj A_sub_mat w_idx w_min A_sub all_edges_rand all_weights_rand edges_removed counter

    end

    WR_tot = [WR_tot min(WR(:,ss))];
    WR_rand_tot = [WR_rand_tot min(WR_rand(:,ss))]; 

end

out_path=('/nfsd/biopetmri4/Users/ClaudiaTarricone/EC_modes_hierarchies/Results/Null_model');
save(fullfile(out_path,'Results_test_null_DAG.mat'),'WR_tot','WR_rand_tot','WR_rand','WR')

%%
f=figure
boxplot([WR_tot',WR_rand_tot'],'Labels',{'original','randomized'})
title('Weights removed (144 subjects)')
set(f, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

path_save=('/nfsd/biopetmri4/Users/ClaudiaTarricone/EC_modes_hierarchies/Results/Null_model');

saveas(f,fullfile(path_save,['Comparison_removed_edges_null_model.png'])) 
%save('...','DAG','MATRIX_EC_ANTISYM','MATRIX_S',...
%    'ER','ER_perc','cycles_rmd')
