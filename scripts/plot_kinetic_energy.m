clear all
clc
close all

addpath('/nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/FINAL_VERSION_CODE/functions')

%% load data (load struct with 0 Eim but figure is the same)
load /nfsd/biopetmri4/Projects/lavori_MRI/HCP_sDCM/Codes/EC_MODES/DATA/S_Sigma_HCP_110625.mat

%%  Kinetic energy of each subject-RUN 1

% Real
for sbj=1:numel(E_tot)-1

    %linear interpolation
    x_interp=[-0.7:0.025:-0.05];
    tmp=griddedInterpolant(E_tot(sbj).gamma_tot_run1,E_tot(sbj).E_tot_run1,'linear','nearest');
    E_re_interp(:,sbj)=tmp(x_interp);

%     figure(1)
%     plot(E_tot(sbj).gamma_tot_run1,E_tot(sbj).E_tot_run1,'r.-')
%     hold on
%     xline(E_tot(sbj).cut_re_run1,Label='energy cutoff')
%     ylabel('Kinetic energy (%)')
%     xlabel('real \lambda')
%     title(['Kinetic energy-sbj=', num2str(sbj)])
%     xlim([-0.7 0])
%     hold off
%     pause()
end

% Imag
for sbj=1:numel(E_tot)-1

    %linear interpolation
    y_interp=[0:0.025:0.4];
    tmp=griddedInterpolant(unique(E_tot(sbj).w_tot_run1),unique(E_tot(sbj).E_tot_run1_Im),'linear','nearest');
    E_im_interp(:,sbj)=tmp(y_interp);

%     figure(1)
%     plot(E_tot(sbj).w_tot_run1,E_tot(sbj).E_tot_run1_Im,'r.-')
%     hold on
%     xline(E_tot(sbj).cut_im_run1,Label='energy cutoff')
%     ylabel('Kinetic energy (%)')
%     xlabel('Im \lambda')
%     title(['Kinetic energy-sbj=', num2str(sbj)])
%     hold off
%     pause()
end

%group-level figure
subplot(211)
stdshade(E_re_interp',0.25,'black',x_interp,0.2)
hold on
errorbar(x_interp,mean(E_re_interp'),std(E_re_interp'),'.-','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
xline(mean([E_tot.cut_re_run1]),Label='energy cutoff')
hold on
xline(mean([E_tot.cut_re_run1])+std([E_tot.cut_re_run1]),'k--')
hold on
xline(mean([E_tot.cut_re_run1])-std([E_tot.cut_re_run1]),'k--')
xlim([min(x_interp) max(x_interp)])
ylabel('Kinetic real energy')
xlabel('real \lambda')
%set patches
% hold on
% patch([min(x_interp) max(x_interp) max(x_interp) min(x_interp)],[-0.1 -0.1 thr1 thr1],[0.8510    0.3255    0.0980],...
%     'FaceAlpha',0.2)
% hold on
% patch([min(x_interp) max(x_interp) max(x_interp) min(x_interp)],[thr1 thr1 thr2 thr2],[0.9294    0.6941    0.1255],...
%     'FaceAlpha',0.2)
% hold on
% patch([min(x_interp) max(x_interp) max(x_interp) min(x_interp)],[thr2 thr2 105 105],[0    0.4471    0.7412],...
%     'FaceAlpha',0.2)
% title('Dissipative kinetic energy')
% hold off

subplot(212)
stdshade(E_im_interp',0.25,'black',y_interp,0.2)
hold on
errorbar(y_interp,mean(E_im_interp'),std(E_im_interp'),'.-','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
xline(mean([E_tot.cut_im_run1]),Label='energy cutoff')
hold on
xline(mean([E_tot.cut_im_run1])+std([E_tot.cut_im_run1]),'k--')
hold on
xline(mean([E_tot.cut_im_run1])-std([E_tot.cut_im_run1]),'k--')
xlim([min(y_interp) max(y_interp)])
ylabel('Kinetic imag energy')
xlabel('imag \lambda')
sgtitle('RUN 1')




