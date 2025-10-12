function CONN_ICC_A = f_ICC_edgewise(Atest,Aretest)
% Evaluate ICC on the connICA matrices
% Add ICC & Network Scripts to path
% addpath(fullfile(pwd,'ICC'))

numEdges = size(Atest,1);
CONN_ICC_A = zeros(1,numEdges);        

for comp=1:numEdges
    data4icc_A = [Atest(comp,:)' Aretest(comp,:)'];
    rows2delete_A = isnan(sum(data4icc_A,2));
    data4icc_A(rows2delete_A,:) = [];
    CONN_ICC_A(1,comp) = ICC(data4icc_A,'1-1') ; 
end

return;