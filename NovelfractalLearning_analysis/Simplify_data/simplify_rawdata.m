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
% %Lemmy
% filepath = 'Y:\PLEXON_GRAYARRAY_LEMMY\NovelfractalLearning_Lemmy_combinedsorting';
% DDD = dir(fullfile(filepath,'*tasksession_mega_file.mat'));
% %  new file location:
% filesavepath = 'Y:\PLEXON_GRAYARRAY_LEMMY\NovelfractalLearning_Lemmy_for_share';


for ii = 1:numel(DDD)
    clear Generaltask SPK* Crossval_neuronlist Neuronlist
    load(fullfile(DDD(ii).folder, DDD(ii).name), 'Generaltask', 'SPK*', 'Crossval_neuronlist');
    
    Neuronlist = Crossval_neuronlist{1,1};
    allfieldnames = fieldnames(Neuronlist);
    fieldname_2_include = {'name';'electrodeID'; 'ISIviolationrate'; 'unitQuality'; 'contaminationRate'; 'exclude'};    
    fieldname_2_remove = setdiff(allfieldnames, fieldname_2_include);
    % Save minimum info in Neuronlist
    Neuronlist([Neuronlist.exclude]==1)=[];
    Neuronlist = rmfield(Neuronlist,fieldname_2_remove);
    
    
    % Save minimum info in Generaltask too
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
    if numel(who('SPK*'))>0
        save(fullfile(filesavepath, DDD(ii).name), 'Generaltask', 'SPK*', 'Neuronlist');
    else
        save(fullfile(filesavepath, DDD(ii).name), 'Generaltask', 'Neuronlist');
    end
end