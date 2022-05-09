%%% This script is aimed for simplifying the data structure in The TSMF,
%%% only keep the necessary data for the plots in the figure and remove all
%%% the rest.
% spike data: ISIviolationrate, unitQuality, contaminationRate, fracsptimes, sptimes,
% electrodeID, electrodeDepth
close all; clear all;

%Slayer original file location:
filepath = 'Y:\PLEXON_GRAYARRAY_Slayer\NovelFractalLearning_Slayer_combinedsorting';
DDD = dir(fullfile(filepath,'*tasksession_mega_file.mat'));
% new file location:
filesavepath = 'Y:\PLEXON_GRAYARRAY_Slayer\NovelFractalLearning_Slayer_for_share';
%Lemmy
% filepath = 'Y:\PLEXON_GRAYARRAY_LEMMY\NovelfractalLearning_Lemmy_combinedsorting';
% DDD = dir(fullfile(filepath,'*tasksession_mega_file.mat'));
% %  new file location:
% filesavepath = 'Y:\PLEXON_GRAYARRAY_LEMMY\NovelfractalLearning_Lemmy_for_share';

% load the neuronlist that combines the infotask result

load(fullfile('Y:\Kaining\NovelfractalLearning_analysis\NovelfractalLearning_analysis\Plot_all_session_code','neuronliststruct_all'), 'Neuronlist_good');
filenamelist = {Neuronlist_good.filename};

%sum of the neurons
n_neurons = 0;

for ii = 1:numel(DDD)
    %profile on
    clear Generaltask SPK* Crossval_neuronlist Neuronlist
    load(fullfile(DDD(ii).folder, DDD(ii).name), 'Generaltask', 'SPK*', 'Crossval_neuronlist');
    
    Neuronlist = Crossval_neuronlist{1,1};
    allfieldnames = fieldnames(Neuronlist);
    fieldname_2_include = {'name';'electrodeID'; 'ISIviolationrate'; 'unitQuality'; 'contaminationRate'; 'exclude';...
        'recency_ind_match_pos'; 'P_recency_ind_match_pos'; 'pred_nov_vs_fam'; 'P_pred_nov_vs_fam'; 'pred_vs_unpred_fam'; 'P_pred_vs_unpred_fam_perm';...
        'All_Fractal_FR'};    
    fieldname_2_remove = setdiff(allfieldnames, fieldname_2_include);
    % Save minimum info in Neuronlist
    Neuronlist = rmfield(Neuronlist,fieldname_2_remove);
    % Remove the excluded neurons
    currentfileInd = strcmp(filenamelist, DDD(ii).name);
    Neuronlist_local = Neuronlist_good(currentfileInd);
    good_neuron = {Neuronlist_local.name};
    good_ind = contains({Neuronlist.name}, good_neuron);
    Neuronlist = Neuronlist(good_ind);
    %Neuronlist([Neuronlist.exclude]==1)=[];
    
    
    % Save minimum info in Generaltask
    fieldname_2_include_Generaltask = {'trialstart'; 'trialend'; 'successtrial'; 'successtrstart'; 'trialtype';...
        'trialnumber'; 'timetargeton'; 'Set'; 'timetargetoff'; 'timefpon'; 'timefpoff'; 'rewardeliver'; 'rewfixangle';...
        'fixacq'; 'b_error_place';  'novelfractaldate'; 'Fractals'};
    allfieldnames_Generaltask = fieldnames(Generaltask);
    fieldname_2_remove_Generaltask = setdiff(allfieldnames_Generaltask, fieldname_2_include_Generaltask);
    Generaltask = rmfield(Generaltask,fieldname_2_remove_Generaltask);
    
    
    % clear the rest of the non-included neurons
    spikenames = who('SPK*');
    excludevarnames = setdiff(spikenames, {Neuronlist.name});
    if ~isempty(excludevarnames)
        clear(excludevarnames{:});
    end
    spikenames = who('SPK*');
    
    % For each spike, only include the necessary neurons
    fieldname_2_include_SPK = {'name';'electrodeID'; 'ISIviolationrate'; 'unitQuality'; 'contaminationRate'; 'sptimes'; 'fracsptimes'};
    
    for iii = 1: numel(spikenames)
        eval(['current_neuron = ' spikenames{iii} ';']);
        allfieldnames_SPK = fieldnames(current_neuron);
        fieldname_2_remove_SPK = setdiff(allfieldnames_SPK, fieldname_2_include_SPK);
        current_neuron = rmfield(current_neuron,fieldname_2_remove_SPK);
        eval([spikenames{iii} '= current_neuron;' ]);
    end
    
    % Check that the Neuronlist matches the SPK to be saved.
    assert(numel(who('SPK*'))==numel(Neuronlist), 'The Neuronlist do not match the raw neurons')
    
    old_filename = DDD(ii).name;
    endind = strfind(old_filename,'_tasksession')-1;
    startind = 10;
    new_filename = ['task_session_', old_filename(startind:endind) '.mat'];
    
    if numel(who('SPK*'))>0
        save(fullfile(filesavepath, new_filename), 'Generaltask', 'SPK*', 'Neuronlist');
    else
        save(fullfile(filesavepath, new_filename), 'Generaltask', 'Neuronlist');
    end
    
    n_neurons = n_neurons+numel(Neuronlist);
    %profile viewer
end

disp(n_neurons);

