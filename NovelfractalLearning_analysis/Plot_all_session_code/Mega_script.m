clear all; clc; close all; restoredefaultpath;
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
load('./maindata/neuronliststruct_Monkey_S.mat');
%Monkey L
load('./maindata/neuronliststruct_Monkey_L.mat');
Neuronlist_all = [Neuronlist_S, Neuronlist_L];
clear Neuronlist_S Neuronlist_L


clear indices
% grab the key indices for hypotheses testing
indices.pred_nov_vs_fam = [Neuronlist_all(:).pred_nov_vs_fam]';
indices.pred_vs_unpred_fam=[Neuronlist_all(:).pred_vs_unpred_fam]';
indices.violation_ind=[Neuronlist_all(:).violation_ind]';
indices.recency_ind=[Neuronlist_all(:).recency_ind_match_pos]';
indices.rewardvalueindex = [Neuronlist_all(:).rewardvalueindex_precue]';
indices.RewInfoAnticipIndex = [Neuronlist_all(:).RewInfoAnticipIndex_split]';

% grab the key indices p values

indices.Ppred_nov_vs_fam = [Neuronlist_all(:).P_pred_nov_vs_fam]';
indices.Ppred_vs_unpred_fam=[Neuronlist_all(:).P_pred_vs_unpred_fam_perm]';
indices.Pviolation_ind=[Neuronlist_all(:).P_violation_ind_perm]';
indices.Precency_ind=[Neuronlist_all(:).P_recency_ind_match_pos]';
indices.Prewardvalueindex = [Neuronlist_all(:).rewardvalueindexP_precue]';
indices.PRewInfoAnticipIndex = [Neuronlist_all(:).RewInfoAnticipIndexP_split]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Analysis start

addpath('./population_analyses');

% Bar plot figure

plotname = ['Indices_barplot_all_session.pdf'];
barplot_func(indices, plotpath, plotname);

% heatmap figure
correlation_plot_func(indices, plotpath, {});

% Info plot
infonew_correlation_barplots;

% control analysis for correlation bar plot
control_analysis_correlation;

%% additional analysis
addpath('./additional_analyses');

mega_additional;

%% Brain area analyses
addpath('./brain_areas_analyses')
load('./brain_areas_analyses/segments/segmentmap')

All_region_plot;

Region_clustering;

barplot_learning_forgetting;


%% Learning analyses
addpath('./learning_analyses');
plotpath = './plots';

Learning_acrossday_2days_timescale;

Learning_acrossday_2days_barplot;

Learning_acrossday_5days_barplot;

fast_slow_learning_forgetting_analysis;


%% pupil
addpath('./pupil');
Monkey = 'L';
load('pupillist_Monkey_L');
pupil_plots(pupillist, Monkey, plotpath);

Monkey = 'S';
load('pupillist_Monkey_S');
pupil_plots(pupillist, Monkey, plotpath);

%% Example neurons:
% Sampleneuron: This requires loading the raw files.
% which can be downloaded here:
% https://wustl.box.com/s/v4x3zjvyopexyud3ghnk87apav9ma6ay
% This code is also an example to show how we calculate novelty index,
% sensory surprise index and recency index.
addpath('./example_neuron_code');
Sampleneuron;

Sampleneuron_for_learning;


profile viewer;

