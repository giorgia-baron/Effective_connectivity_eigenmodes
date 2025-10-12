function num_couple_eig= couple_eig(eigval_sub,N_rois)
    %evaluate couples of eigs
    num_couple_eig=[];
    counter=1;

    while counter<N_rois

        if and(real(eigval_sub(counter))==real(eigval_sub(counter+1)),...
                abs(imag(eigval_sub(counter)))==abs(imag(eigval_sub(counter+1))))
            num_couple_eig=[num_couple_eig 2];
            counter=counter+2;
        else
            num_couple_eig=[num_couple_eig 1];
            counter=counter+1;
        end
            
    end

    if not(sum(num_couple_eig)==N_rois)
        num_couple_eig(end+1)=1;
    end
end