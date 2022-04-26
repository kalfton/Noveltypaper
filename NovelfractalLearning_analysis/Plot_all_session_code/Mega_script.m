clear all; clc; warning off; beep off; close all;
addpath('./help_func');
addpath('./utils');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile on;
plotpath = './plots';
exclude_mua = true;
shuffling_num = 1000; % the shuffling number is set to 1000 to reduce the running time, it is originally set as 10000
StatisticalThreshold=0.01;

% load data:
%Monkey S
load('./maindata/neuronliststruct_Slayer_simplified.mat');
%Monkey L
load('./maindata/neuronliststruct_Lemmy_simplified.mat');
Neuronlist_good = [Neuronlist_Slayer, Neuronlist_Lemmy];
clear Neuronlist_Slayer Neuronlist_Lemmy


clear indices
% grab the key indices for hypotheses testing
indices.pred_nov_vs_fam = [Neuronlist_good(:).pred_nov_vs_fam]';
indices.pred_vs_unpred_fam=[Neuronlist_good(:).pred_vs_unpred_fam]';
indices.violation_ind=[Neuronlist_good(:).violation_ind]';
indices.recency_ind=[Neuronlist_good(:).recency_ind_match_pos]';
indices.rewardvalueindex = [Neuronlist_good(:).rewardvalueindex_precue]';
indices.RewInfoAnticipIndex = [Neuronlist_good(:).RewInfoAnticipIndex_split]';

% grab the key indices p values

indices.Ppred_nov_vs_fam = [Neuronlist_good(:).P_pred_nov_vs_fam]';
indices.Ppred_vs_unpred_fam=[Neuronlist_good(:).P_pred_vs_unpred_fam_perm]';
indices.Pviolation_ind=[Neuronlist_good(:).P_violation_ind_perm]';
indices.Precency_ind=[Neuronlist_good(:).P_recency_ind_match_pos]';
indices.Prewardvalueindex = [Neuronlist_good(:).rewardvalueindexP_precue]';
indices.PRewInfoAnticipIndex = [Neuronlist_good(:).RewInfoAnticipIndexP_split]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% odd trials
% indices.pred_nov_vs_fam_odd = [Neuronlist_odd(:).pred_nov_vs_fam]';
% indices.pred_vs_unpred_fam_odd=[Neuronlist_odd(:).pred_vs_unpred_fam]';
% indices.violation_ind_odd=[Neuronlist_odd(:).violation_ind]';
% indices.recency_ind_odd=[Neuronlist_odd(:).recency_ind_match_pos]';
% 
% indices.Ppred_nov_vs_fam_odd = [Neuronlist_odd(:).P_pred_nov_vs_fam]';
% indices.Ppred_vs_unpred_fam_odd=[Neuronlist_odd(:).P_pred_vs_unpred_fam_perm]';
% indices.Pviolation_ind_odd=[Neuronlist_odd(:).P_violation_ind_perm]';
% indices.Precency_ind_odd=[Neuronlist_odd(:).P_recency_ind_match_pos]';
% 
% 
% %% even trials
% indices.pred_nov_vs_fam_even = [Neuronlist_even(:).pred_nov_vs_fam]';
% indices.pred_vs_unpred_fam_even=[Neuronlist_even(:).pred_vs_unpred_fam]';
% indices.violation_ind_even=[Neuronlist_even(:).violation_ind]';
% indices.recency_ind_even=[Neuronlist_even(:).recency_ind_match_pos]';
% 
% indices.Ppred_nov_vs_fam_even = [Neuronlist_even(:).P_pred_nov_vs_fam]';
% indices.Ppred_vs_unpred_fam_even=[Neuronlist_even(:).P_pred_vs_unpred_fam_perm]';
% indices.Pviolation_ind_even=[Neuronlist_even(:).P_violation_ind_perm]';
% indices.Precency_ind_even=[Neuronlist_even(:).P_recency_ind_match_pos]';

%%

addpath('./population_analyses');

% Bar plot figure

plotname = ['Indices_barplot_all_session.pdf'];
barplot_func(indices, plotpath, plotname);

% heatmap figure
heatmapplot_func(indices, plotpath, {}, shuffling_num);

% SDF figure
tic
if isfield(Neuronlist_good, 'pred_nov_vs_fam_sdfs')
    %SDF_plot_Crossvalidation_V2
    %SDF_plot_Crossvalidation_infonewtask
end
toc

% Info plot
infonew_correlation_barplots;

% 3D space of Novelty neurons
Neuron_3D_space_scatter;

% control analysis for correlation bar plot
control_analysis_correlation;


% example neurons: This need the SDFs.
%Sampleneuron;

%% additional analysis
addpath('./additional_analyses');

mega_additional;

%% Brain area analyses
addpath('./brain_areas_analyses')
load('./brain_areas_analyses/segments/segmentmap')
All_region_plot;

Region_clustering;

barplot_learning_forgetting;

% this is an additional fun analysis related to Murray et. al., 2014
time_scale_check;


%% Learning analyses
addpath('./learning_analyses');
plotpath = './plots';

Learning_acrossday_2days_timescale;

Learning_acrossday_2days_barplot;

% 5 days learning plots take a long time to run
Learning_acrossday_5days_barplot;

fast_slow_learning_forgetting_analysis;

% sample neuron plot need the task information from origin files
% Sampleneuron_for_learning;


%% pupil
addpath('./pupil');
Monkey = 'Lemmy';
load('pupillist_Lemmy');
barplot_pupil_SDF;

Monkey = 'Slayer';
load('pupillist_Slayer');
barplot_pupil_SDF;

profile viewer;

