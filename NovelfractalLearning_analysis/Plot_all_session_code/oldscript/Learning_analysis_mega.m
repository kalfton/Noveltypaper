% Learning analysis mega
addpath('./util_learning');
plotpath = './plots';

%load('neuronliststruct');

Learning_acrossday_2days_timescale;

Learning_acrossday_2days_barplot;

% 5 days learning plots take a long time to run
Learning_acrossday_5days_barplot;

fast_slow_learning_forgetting_analysis;