function [A_sparse,mask_sparse] = sparsify(A,num_rois)
   
%fit bimodal gaussians
mask_A=ones(num_rois)-eye(num_rois);
mask_A=triu(mask_A,1);
A_abs=abs(A-diag(diag(A)));
A_hist=log(A_abs(logical(mask_A)));
[b, x_bin] = histcounts(A_hist);
x_bin = x_bin(1:end-1); %for plotting
%     figure()
%     plot(x_bin, b)
%     hold on

AIC = zeros(1,2);
gmm_model = cell(1,2);
for k = 1:2
    gmm_model{k} = fitgmdist(A_hist, k);
    AIC(k) = gmm_model{k}.AIC;
end
[minAIC, numComp] = min(AIC);
gmm_best = gmm_model{numComp};
if numComp == 1
    area_overlap= [1, 1];
else
    %x_bin = -15:.1:0;
    pdf1 = normpdf(x_bin, gmm_best.mu(1), gmm_best.Sigma(1,1,1)) * gmm_best.ComponentProportion(1); %componentproportion=percentage of area of gaussian to be fitted
    pdf2 = normpdf(x_bin, gmm_best.mu(2), gmm_best.Sigma(1,1,2)) * gmm_best.ComponentProportion(2);
    [~, idx_mu] = sort(gmm_best.mu);
    area_overlap_i = sum(min(pdf1, pdf2)) * (x_bin(2) - x_bin(1)) ./ gmm_best.ComponentProportion(idx_mu);
    area_overlap= area_overlap_i;
    plot(x_bin, pdf1)
    hold on
    plot(x_bin, pdf2)
    hold on
    plot(x_bin, gmm_best.pdf(x_bin'))
    title(['w overlap ', num2str(area_overlap_i)])
    hold on

    %find minimum of distribution
    [max_x1,idx]=max(pdf1);
    x1=find(x_bin==x_bin(idx));
    [max_x2,idx]=max(pdf2);
    x2=find(x_bin==x_bin(idx));
    pdf_biv=gmm_best.pdf(x_bin');
    if x1<x2
        [min_val,min_idx]=min(pdf_biv(x1:x2));
    else
        [min_val,min_idx]=min(pdf_biv(x2:x1));
    end
%         if min_val>gmm_best.mu(2) && min_val<gmm_best.mu(1)
    min_thr=x_bin(pdf_biv==min_val);
    thr=exp(min_thr);
    plot(min_thr,min_val,'ro')
    hold off
    mask_sparse=abs(A)<thr;
    A_sparse=A;
    A_sparse(mask_sparse)=0;
%         end
end
end