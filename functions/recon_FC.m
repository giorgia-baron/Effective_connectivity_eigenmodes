function FC_recon= recon_FC(V,eigval,Q_noiseVar,range)
    %compute FC of each mode: the computation should be equal to that
    %proposed by Giacomo
%     W=inv(V');
%     W=W';
%     w=W(1,:);
%     V_inv_T=(inv(V))';
%     M=W*V_inv_T;
%     M=Q_noiseVar(1,1)*M; -->M should equals Q

            
    T_mm=V(:,range);
    T_mm_inv_all=inv(V);
    T_mm_inv=T_mm_inv_all(range,:);
    T_eig_mm=diag(eigval(range));

    %apply lyap eq
    A=T_eig_mm;
    Q=T_mm_inv*Q_noiseVar*T_mm_inv';
    Sigma_mm=lyap(A,Q);
    FC_recon=Q_noiseVar(1,1)*T_mm*Sigma_mm*T_mm';%Q_noiseVar(1,1)*T_mm*Sigma_mm*T_mm';
    Q_recon=T_mm*Q*T_mm';
end
