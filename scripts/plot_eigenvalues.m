clearvars -except Div_ranges
clc
close all


%% load data (load struct with 0 Eim but figure is the same)
load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/S_Sigma_HCP_110625.mat

%% run1 eigenvalues plot

% Color map 
colorMap = [0    0.4471    0.7412; 0 1 0; 0.9294    0.6941    0.1255; 0.8510    0.3255    0.0980];

sbj_run1=numel(mask_eig_all_sub);

for sbj=1:sbj_run1-1

    clear leng div
    
    %check if real and imag are correctly paired
    if sum(mask_eig_all_sub(sbj).mask_eig_run1)~=2
        disp('error range inconsistency')
        pause()
    end
    
    %length=all complex eigs >-0.7
    Re_low_sbj=mask_eig_all_sub(sbj).mask_eig_run1(1,:);
    Re_high_sbj=mask_eig_all_sub(sbj).mask_eig_run1(2,:);
    Im_low_sbj=mask_eig_all_sub(sbj).mask_eig_run1(3,:);
    Im_high_sbj=mask_eig_all_sub(sbj).mask_eig_run1(4,:);

    for leng=1:length(Re_low_sbj)

        % create 4 ranges from Im/Re combination
        if Re_low_sbj(leng)==1 && Im_low_sbj(leng)==1
            div(leng)=1;
        end
        if Re_low_sbj(leng)==1 && Im_high_sbj(leng)==1
            div(leng)=2;
        end
        if Re_high_sbj(leng)==1 && Im_low_sbj(leng)==1
            div(leng)=3;
        end
        if Re_high_sbj(leng)==1 && Im_high_sbj(leng)==1
            div(leng)=4;
        end
    end

    Div_ranges(sbj).div_range_run1=div;

end


figure;
hold on;

for sbj = 1:sbj_run1-1
    eigs = eigval_all_sub(sbj).eigenall_run1(:); %all complex eigs >-0.7
    
    % eigenvalue classification
    for jj = 1:length(eigs)
        z = eigs(jj);
        
        if not(length(Div_ranges(sbj).div_range_run1)==length(eigs))
            disp('error')
            pause()
        end
        % color assignment
        if Div_ranges(sbj).div_range_run1(jj)==1
            color = colorMap(1, :); % Range 1: blue
        elseif Div_ranges(sbj).div_range_run1(jj)==2
            color = colorMap(2, :); % Range 2: green
        elseif Div_ranges(sbj).div_range_run1(jj)==3
            color = colorMap(3, :); % Range 3: yellow
        elseif Div_ranges(sbj).div_range_run1(jj)==4
            color = colorMap(4, :); % Range 4: red
        end
        
        % eigenvalue plot
        plot(real(z), imag(z), 'o', 'Color', color,'MarkerFaceColor',color,'MarkerSize',5,'MarkerEdgeColor','k');
        xlim([-0.7 0])
        hold on
%         pause();
    end
end

hold off;
xlabel('real(\lambda)');
ylabel('imag(\lambda)');
title('EC eigenvalues');
subtitle('Dividing Re and Im energy components-run1')
grid on;

% %% run1 eigenvalues plot
% 
% % Color map 
% colorMap = [0    0.4471    0.7412; 0 1 0; 0.9294    0.6941    0.1255; 0.8510    0.3255    0.0980];
% 
% sbj_run1=numel(mask_eig_all_sub);
% 
% for sbj=1:sbj_run1-1
% 
%     clear leng div
% 
%     %check if real and imag are correctly paired
%     if sum(mask_eig_all_sub(sbj).mask_eig_run1)~=2
%         disp('error range inconsistency')
%         pause
%     end
%     
%     %length=all complex eigs >-0.7
%     Re_low_sbj=mask_eig_all_sub(sbj).mask_eig_run1(1,:);
%     Re_high_sbj=mask_eig_all_sub(sbj).mask_eig_run1(2,:);
%     Im_low_sbj=mask_eig_all_sub(sbj).mask_eig_run1(3,:);
%     Im_high_sbj=mask_eig_all_sub(sbj).mask_eig_run1(4,:);
% 
%     for leng=1:length(Re_low_sbj)
% 
%         % create 4 ranges from Im/Re combination
% 
%         if Re_low_sbj(leng)==1 && Im_low_sbj(leng)==1
%             div(leng)=1;
%         end
%         if Re_low_sbj(leng)==1 && Im_high_sbj(leng)==1
%             div(leng)=2;
%         end
%         if Re_high_sbj(leng)==1 && Im_low_sbj(leng)==1
%             div(leng)=3;
%         end
%         if Re_high_sbj(leng)==1 && Im_high_sbj(leng)==1
%             div(leng)=4;
%         end
%     end
% 
% Div_ranges(sbj).div_range_run1=div;
% 
% end
% 
% 
% figure;
% hold on;
% 
% for sbj = 1:sbj_run1-1
%     eigs = eigval_all_sub(sbj).eigenall_run1(:); %all complex eigs >-0.7
%     
%     % eigenvalue classification
%     for jj = 1:length(eigs)
%         z = eigs(jj);
%         
%         if not(length(Div_ranges(sbj).div_range_run1)==length(eigs))
%             disp('error')
%             pause
%         end
%         % Color assignment
%         if Div_ranges(sbj).div_range_run1(jj)==1
%             color = colorMap(1, :); % Range 1: colore blu
%         elseif Div_ranges(sbj).div_range_run1(jj)==2
%             color = colorMap(2, :); % Range 2: colore verde
%         elseif Div_ranges(sbj).div_range_run1(jj)==3
%             color = colorMap(3, :); % Range 3: colore arancione
%         elseif Div_ranges(sbj).div_range_run1(jj)==4
%             color = colorMap(4, :); % Range 4: colore rosso
%         end
%         
%         % eigenvalue plot
%         plot(real(z), imag(z), 'o', 'Color', color,'MarkerFaceColor',color,'MarkerSize',5,'MarkerEdgeColor','k');
% %         plot(real(z), -imag(z), 'o', 'Color', color,'MarkerFaceColor',color,'MarkerSize',5); % Simmetria per i coniugati
%         xlim([-0.7 0])
%         hold on
% %         pause();
%     end
% end
% 
% hold off;
% xlabel('real(\lambda)');
% ylabel('imag(\lambda)');
% title('EC eigenvalues');
% subtitle('Dividing Re and Im energy components-run1')
% grid on;
% 
% 
