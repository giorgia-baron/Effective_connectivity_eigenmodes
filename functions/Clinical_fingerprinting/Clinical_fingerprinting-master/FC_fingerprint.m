%% Batch example of the Clinical Connectome Fingerprint analysis, proposed in  (Sorrentino P. et al., Neuroimage 2021),
% to compute identifiability in clinical populations.
% The code compares two functional connectome (FC) acquisitions/sessions for each subject (defined as test and retest)
% contained in the sample dataset (data_test_FC) and give as result an
% "identifiabilty matrix" (Amico & Goñi, Scientific Reports 2018). Its main diagonal (top left to bottom right) highlights the
% self identifiability (Iself) of the sample, i.e. the similarity between the test and retest FC of the
% same subject. The remaining elements consists in the comparison between
% test and retest FC of different subjects (Iothers).
% Clinical connetome fingerprinting extends this idea by comparing
% individual connectomes between different populations (i.e. Controls and Mild Cognitive
% Impairment, see Sorrentino et al. for details).
% Functional Connectomes matrices are obtained through Phase Linearity
% Measurement (PLM) (Baselice et al., 2019), here only reported in the Alpha band.

%% PLEASE NOTE: For privacy reasons, we could not use the same data as in Sorrentino et al. 
%% Hence, the sample results reported here do not correspond to the original manuscript, instead they are merely illustrative of the methodology.
%% Please contact giuseppe.sorrentino@uniparthenope.it for 
%% (reasonably) requesting the original sample connectomes of the manuscript.
%% The data used in this code come from healthy MEG connectomes obtained from Human Connectome Project dataset (see Sareen et al., NeuroImage 2021, for details).

%% Authors: Emahnuel TROISI LOPEZ, Pierpaolo SORRENTINO & Enrico AMICO 
% version 1.0. July 01, 2021
%
%% PLEASE CITE US!
% If you are using this code for your research, please kindly cite us:

% Authors: Sorrentino P, Rucco R, Lardone A, Liparoti M, Lopez ET, Cavaliere C, Soricelli A, Jirsa V, Sorrentino G, Amico E.
% Title: Clinical connectome fingerprints of cognitive decline.
% Published on: NeuroImage - 2021 Jun 9, p. 118253.
% Doi: doi.org/10.1016/j.neuroimage.2021.118253

%% Initialize environment
clearvars
clc
close all;


%% Load FC (functional connectivity) resting-state data. This is a sample data.
% loaded data consists in a cell for each subject. Each cell contains a 3D
% matrix: dimensions are brain regions x brain regions x test/re-test FC
load('data_test_FC.mat');

n_subj = size(data_test_FC,2); %Number of subjects
n_roi = size(data_test_FC{1},1); %Number of region of interest (ROI)

%% Start

%Convert each FC matrix into an array. See details in the f_load_mat function
mask_ut = triu(true(n_roi,n_roi),1);
[FCs_test,FCs_retest] = f_load_mat(data_test_FC,mask_ut);
mask_diag = logical(eye(size(FCs_test,1)));

%Identifiability matrix
Ident_mat = corr(FCs_test',FCs_retest');
Ident_mat_1 = Ident_mat;

%Success rate: how many times an Iself value is higher than an Iothers
%value on the same row and column.
sr1 = zeros(1,n_subj);
for i=1:length(Ident_mat)
    sr1(i) = sum(Ident_mat(i,i)>Ident_mat(i,:))+sum(Ident_mat(i,i)>Ident_mat(:,i));
end
sr_indiv = sr1*100/(length(Ident_mat)*2-2); %individual success rate
sr_mean = mean(sr_indiv,2); %average success rate

%Parameters
Iself = zeros(1,n_subj); %self identifiability
Iothers = zeros(1,n_subj); % Others-identifiability
Idiff_indiv = zeros(1,n_subj); %Iself - Iothers
Idiff = nanmean(Ident_mat(mask_diag))-nanmean(Ident_mat(~mask_diag));
for s=1:n_subj
    Iself(s) = Ident_mat(s,s);
    Ident_mat(s,s) = nan;
    Iothers(s) = 0.5.*(nanmean(Ident_mat(s,:)) + nanmean(Ident_mat(:,s))');
    Idiff_indiv(s) = Iself(s)-Iothers(s);
end

%Plot
IDMAT = figure;
subplot(1,2,1)
imagesc(Ident_mat_1); axis square; set(gca,'Xtick',[]);set(gca,'Ytick',[]); colorbar; caxis([0.3 .9]);
title('Identifiability Matrix');
xlabel('Test FC')
ylabel('Retest FC')

%% Perform intraclass correlation to highlight which roi mainly contributes to the identifiability
addpath(fullfile(pwd,'ICC'))
ICC_threshold = 0.68;
disp('Computing ICC..')
ICC_struct = f_ICC_edgewise(FCs_test(:,:)',FCs_retest(:,:)');
ICC_mat = zeros(n_roi,n_roi);
ICC_mat(mask_ut) = ICC_struct;
ICC_mat = ICC_mat + ICC_mat';

subplot(1,2,2)
imagesc(ICC_mat); axis square; set(gca,'Xtick',[]);set(gca,'Ytick',[]); colorbar; %caxis([-.1 .8]);
title('ICC Matrix');
xlabel('90 roi')
ylabel('90 roi')

%% Iclinical
%Here we illustrate how the Iclinical score is computed and how an
%Iclinical matrix appears. For this purpose we build two sample groups
%from our dataset. Now we will create 10 test and retest FCs for controls and 10
%test and retest FCs to simulate a clinical group.

%Control group
sub_FCs_test_ctrl = FCs_test(1:10,:);
sub_FCs_retest_ctrl = FCs_retest(1:10,:);

%Simulated clinical group
sub_FCs_test_clin = FCs_test(11:20,:);
sub_FCs_retest_clin = FCs_retest(11:20,:);

%Iclinical matrices are obtained cross correlating test of controls with
%retest of clinicals and vice versa.
Iclinical_mat1 = corr(sub_FCs_test_ctrl',sub_FCs_retest_clin');
Iclinical_mat2 = corr(sub_FCs_test_clin',sub_FCs_retest_ctrl');

%Finally, the Iclinical score is computed as the mean similarity of each
%clinical subject with each control subject, averaged between the two
%Iclinical matrices.
Iclinical1=nanmean(Iclinical_mat1,1);
Iclinical2=nanmean(Iclinical_mat2,2)';
Iclinical=nanmean([Iclinical1;Iclinical2],1);

ICMAT=figure;
subplot(1,2,1)
imagesc(Iclinical_mat1); axis square; set(gca,'Xtick',[]);set(gca,'Ytick',[]); colorbar; caxis([0.3 .9]);
xlabel('Controls´ FC test')
ylabel('Clinicals´ FC retest')
title('Iclinical Matrix #1')

subplot(1,2,2)
imagesc(Iclinical_mat2); axis square; set(gca,'Xtick',[]);set(gca,'Ytick',[]); colorbar; caxis([0.3 .9]);
xlabel('Clinicals´ FC test')
ylabel('Controls´ FC retest')
title('Iclinical Matrix #2')


