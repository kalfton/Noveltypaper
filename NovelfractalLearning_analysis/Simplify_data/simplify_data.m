%%% This script is aimed for simplifying the data structure in The TSMF,
%%% only keep the necessary data for the plots in the figure and remove all
%%% the rest.
% spike data: ISIviolationrate, unitQuality, contaminationRate, fracsptimes, sptimes,
% electrodeID, electrodeDepth

% LFP: none; MUA none;

% AI01-05: yes (what is AI05? remove if unnecessary. )

% Generaltask: delete the names of the other tasks and AI and LFP.
% Crossval_neuronlist: (without the _new!!!) Delete pupil data(?)
% Neuronlist in the coressval_neuronlist: remove MUA and exclude neurons!
% Maybe I don't even need to save the crossval, just save the original
% neuronlist wich include:
% name electrodeID electrodeDepth ISIviolationrate, unitQuality, contaminationRate
% 



% Save the Neuronlist/Neuronlist_odd/Neuronlist_even of all recordings, such that the code for plots don't need to
% rerun the indices.
% fields to save: electrodeID, electrodeDepth, ISIviolationrate, unitQuality, contaminationRate
% recency_ind_match_pos, P_recency_ind_match_pos, pred_nov_vs_fam, P_pred_nov_vs_fam,
% pred_vs_unpred_fam, P_pred_vs_unpred_fam_perm, violation_ind, P_violation_ind_perm
% pred_nov_vs_fam_control1, pred_nov_vs_fam_control2,
% P_pred_nov_vs_fam_control1, P_pred_nov_vs_fam_control2,
% pred_vs_unpred_fam_old, P_pred_vs_unpred_fam_perm_old,
% pred_vs_unpred_fam_control, P_pred_vs_unpred_fam_perm_control,
% filename, All_Fractal_FR, frac_sets, zscorestd, zscoremean, monkeyName,
% MonkeyID, region, regionIndex, 
% rewardvalueindex_precue, rewardvalueindexP_precue,
% RewInfoAnticipIndex_split, RewInfoAnticipIndexP_split
% learning, nonlearning, 'Learning_Familiar_start', 'Learning_Novel_start', 'Learning_1day_start', 'Learning_2day_start'
% 'Learning_Familiar_end', 'Learning_Novel_end', 'Learning_1day_end', 'Learning_2day_end'
% 'Learning_Familiar_startroc', 'Learning_Novel_startroc', 'Learning_1day_startroc', 'Learning_2day_startroc'
% 'Learning_Familiar_endroc', 'Learning_Novel_endroc', 'Learning_1day_endroc', 'Learning_2day_endroc'
% 'Learning_Familiar_startp', 'Learning_Novel_startp', 'Learning_1day_startp', 'Learning_2day_startp'
% 'Learning_Familiar_endp', 'Learning_Novel_endp', 'Learning_1day_endp', 'Learning_2day_endp'
% learningforgetinganalysis, withindaylearning_newindex, acrossdayforget_newindex
% withindaylearning_newindex_p, acrossdayforget_newindex_p
% P_object_select_novel_start, obj_select_ind_novel_start;
% P_object_select_novel_end, obj_select_ind_novel_end


% Neuronlist simplification
load(fullfile('Y:\Kaining\NovelfractalLearning_analysis\NovelfractalLearning_analysis\Plot_all_session_code','neuronliststruct_all'), 'Neuronlist_good'); %, 'Neuronlist_even', 'Neuronlist_odd'
allfieldnames = fieldnames(Neuronlist_good);
fieldname_2_include = {'name';'electrodeID'; 'electrodeDepth'; 'fcsvCoordinates_Matlab'; 'ISIviolationrate'; 'unitQuality'; 'contaminationRate'; 'exclude'; 'excludemua';
'recency_ind_match_pos'; 'P_recency_ind_match_pos'; 'pred_nov_vs_fam'; 'P_pred_nov_vs_fam';
'pred_vs_unpred_fam'; 'P_pred_vs_unpred_fam_perm'; 'violation_ind'; 'P_violation_ind_perm';
'pred_nov_vs_fam_control1'; 'pred_nov_vs_fam_control2';
'P_pred_nov_vs_fam_control1'; 'P_pred_nov_vs_fam_control2';
'pred_vs_unpred_fam_old'; 'P_pred_vs_unpred_fam_perm_old';
'pred_vs_unpred_fam_control'; 'P_pred_vs_unpred_fam_perm_control';
'filename'; 'All_Fractal_FR'; 'frac_sets'; 'zscorestd'; 'zscoremean'; 'monkeyName';
'MonkeyID'; 'region'; 'regionIndex'; 
'rewardvalueindex_precue'; 'rewardvalueindexP_precue';
'RewInfoAnticipIndex_split'; 'RewInfoAnticipIndexP_split'
'learning'; 'nonlearning'; 'Learning_Familiar_start'; 'Learning_Novel_start'; 'Learning_1day_start'; 'Learning_2day_start';
'Learning_Familiar_end'; 'Learning_Novel_end'; 'Learning_1day_end'; 'Learning_2day_end';
'Learning_Familiar_startroc'; 'Learning_Novel_startroc'; 'Learning_1day_startroc'; 'Learning_2day_startroc';
'Learning_Familiar_endroc'; 'Learning_Novel_endroc'; 'Learning_1day_endroc'; 'Learning_2day_endroc';
'Learning_Familiar_startp'; 'Learning_Novel_startp'; 'Learning_1day_startp'; 'Learning_2day_startp';
'Learning_Familiar_endp'; 'Learning_Novel_endp'; 'Learning_1day_endp'; 'Learning_2day_endp'; 'learning_surprise';
'learningforgetinganalysis'; 'withindaylearning_newindex'; 'acrossdayforget_newindex';
'withindaylearning_newindex_p'; 'acrossdayforget_newindex_p';
'P_object_select_novel_start'; 'obj_select_ind_novel_start';
'P_object_select_novel_end'; 'obj_select_ind_novel_end';
};

% 'pred_nov_vs_fam_sdfs'; 'recency_ind__match_pos_sdfs'; 'pred_vs_unpred_fam_sdfs'; 'violation_ind_comp_sdf';
% 'SDF_rew100ni'; 'SDF_rew50ni_del'; 'SDF_rew50ni_ndel'; 'SDF_rew0ni';
% 'SDF_rew100i'; 'SDF_rew50i_del'; 'SDF_rew50i_ndel'; 'SDF_rew0i';
% 'pred_nov_vs_fam_sdfs';'recency_ind_sdfs'; 'pred_vs_unpred_fam_sdfs'; 'violation_ind_comp_sdf';
% 'zscoremean_info'; 'zscorestd_info'; 'zscoremean'; 'zscorestd'; 'intrinsic_time';



fieldname_2_remove = setdiff(allfieldnames, fieldname_2_include);

% remove mua and excluded neurons
Neuronlist_good([Neuronlist_good.exclude]==1)=[];
% Neuronlist_even([Neuronlist_good.exclude]==1)=[];
% Neuronlist_odd([Neuronlist_good.exclude]==1)=[];


Neuronlist_good = rmfield(Neuronlist_good,fieldname_2_remove);
% Neuronlist_even = rmfield(Neuronlist_even,fieldname_2_remove(isfield(Neuronlist_even,fieldname_2_remove)));
% Neuronlist_odd = rmfield(Neuronlist_odd,fieldname_2_remove(isfield(Neuronlist_odd,fieldname_2_remove)));


% if the depth is negative or too large make them nan

for ii = 1: numel(Neuronlist_good)
    if Neuronlist_good(ii).electrodeDepth<0 || Neuronlist_good(ii).electrodeDepth>100
        %Neuronlist_good(ii).electrodeDepth = nan;
        Neuronlist_good(ii).fcsvCoordinates_Matlab = [nan, nan, nan];
    end
end
Neuronlist_good = rmfield(Neuronlist_good,'electrodeDepth');


% fix one specific neuron:
neuronind = find(cellfun(@isempty,{Neuronlist_good(:).region}));%8142
%Neuronlist_good(ii).electrodeDepth = nan;
Neuronlist_good(neuronind).fcsvCoordinates_Matlab = [nan, nan, nan];
Neuronlist_good(neuronind).region = "Not Assigned";
Neuronlist_good(neuronind).regionIndex = 45;

Neuronlist_good = RenameField(Neuronlist_good, 'fcsvCoordinates_Matlab', 'electrodeLocation');


% separate Slayer and Lemmy's neuronlist

Neuronlist_Slayer = Neuronlist_good([Neuronlist_good.MonkeyID]==1);
Neuronlist_Lemmy = Neuronlist_good([Neuronlist_good.MonkeyID]==2);
save(fullfile('..\Plot_all_session_code\maindata','neuronliststruct_Lemmy_simplified'), 'Neuronlist_Slayer');
save(fullfile('..\Plot_all_session_code\maindata','neuronliststruct_Slayer_simplified'), 'Neuronlist_Lemmy');

%% 
cd ../Plot_all_session_code
Mega_script;


