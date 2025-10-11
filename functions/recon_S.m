function S_recon = recon_S(V,eigval,Q_noiseVar,range,EC_sub)
    %reconstruct S matrix (complex eigs together)
    %S corresponding to real eigs should equal 0
    
    %V=column right eigenvectors
    %eigval=vector of eigenvalues
    T_mm=V(:,range); %V1 %range determines which mode I'm considering (length(range)=1 or =2)
    T_mm_inv_all=inv(V);
    T_mm_inv=T_mm_inv_all(range,:);
    T_eig_mm=diag(eigval(range)); %L1

    %apply lyap eq
    A=T_eig_mm;
    Q=T_mm_inv*Q_noiseVar*T_mm_inv';
    Sigma_mm1=lyap(A,Q); %Sigmax1
    
    %alternative way to compute Sigmax1
    Sigma_sub=lyap(EC_sub,Q_noiseVar);
    Sigma_mm=T_mm_inv_all*Sigma_sub*T_mm_inv_all';
    Sigma_mm=Sigma_mm(range,range); %checked that it is equal to Sigma_mm1

    S_recon=0.5*(T_mm*T_eig_mm*Sigma_mm1*T_mm'-T_mm*Sigma_mm1*T_eig_mm'*T_mm');  

      

end