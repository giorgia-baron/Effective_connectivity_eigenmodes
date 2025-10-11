function [FCs_test,FCs_retest] = f_load_mat(FCtype,mask_ut)
% Convert matrices into array.
% FCtype is a cell array. Each cell contains two functional connectome (FC)
% of the same subject (name test and retest), according to the following
% dimensions: brain regions x brain regions x test/retest.
% mask_ut is a logical mask build upon the FCs, which only selects the
% upper triangular element of the FC matrices contained in FCtype.
% FCs_test is a matrix composed by a row for each subject and a column for
% each element of the upper triangular of the test FC contained in FCtype.
% FCs_retest is the same as FCs_test, but composed by the elements of the
% retest FC contained in FCtype.

FCs_test = nan(length(FCtype),nnz(mask_ut));
FCs_retest = nan(length(FCtype),nnz(mask_ut));

for i=1:length(FCtype)
    
    aux_test = FCtype{i}(:,:,1);
    aux_retest = FCtype{i}(:,:,2);
    
    FCs_test(i,:) = aux_test(mask_ut);
    FCs_retest(i,:) = aux_retest(mask_ut);
end
