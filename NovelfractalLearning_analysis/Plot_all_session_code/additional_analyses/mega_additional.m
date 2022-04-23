% An infra structure script which integrates all the analysis for reviewer.
addpath('..\utils');
addpath('..\help_func');
rng(0);
DDD=[];
datapath = 'Y:\PLEXON_GRAYARRAY_Slayer\NovelFractalLearning_Slayer_combinedsorting';
DDD = dir(fullfile(datapath,'*tasksession_mega_file.mat'));
datapath = 'Y:\PLEXON_GRAYARRAY_LEMMY\NovelfractalLearning_Lemmy_combinedsorting';
DDD = [DDD; dir(fullfile(datapath,'*tasksession_mega_file.mat'))];

%plotpath = '.\plots';

generate_new_data = false; %true or false
tic
if generate_new_data
    fixation_obj = fixation_session();
    classifer_obj = novelty_classifier_class();
    noise_corr_obj = noise_correlation_class('noveltyresponsive');
    
    
    for fileID = 1:length(DDD)
        clear Generaltask Crossval_neuronlist
        filename = DDD(fileID).name;
        datapath = DDD(fileID).folder;
        
        load(fullfile(datapath, filename),'Generaltask');
        
        load(fullfile(datapath, filename),'Crossval_neuronlist');
        Neuronlist = Crossval_neuronlist{1,1};
        
        % fixation anaylsis here
        clear AI04;
        load(fullfile(datapath, filename),'AI04');
        pupilchannel = AI04;
        
        if contains(datapath, 'Lemmy', 'IgnoreCase', true)
            Monkey = 'Lemmy';
        elseif contains(datapath, 'Slayer', 'IgnoreCase', true)
            Monkey = 'Slayer';
        end
        
        % Check whether there are valid neurons in the session
        if nnz(contains({Neuronlist.name}, 'SPK') & [Neuronlist.exclude]==0)<1
            continue;
        end
        
        fixation_obj.fixation_data_collect(Generaltask);
        
        % classifer analysis here
        classifer_obj.novelty_classifier_singlesession(Generaltask, Neuronlist,0);
        
        % noise correlation analysis here
        noise_corr_obj.noise_correlation_singlesession(Generaltask, Neuronlist);
        
        
        fprintf('%d files done\n', fileID);
    end
    
    % save the data for later quicker use
    propertyname = properties(fixation_obj);
    fixation_struct = struct(); %use a struct to store the information
    for ii = 1:numel(propertyname)
        fixation_struct.(propertyname{ii}) = fixation_obj.(propertyname{ii});
    end
    save('./additional_analyses/data/fixation_struct.mat', 'fixation_struct');
    
    propertyname = properties(classifer_obj);
    classifer_struct = struct(); %use a struct to store the information
    for ii = 1:numel(propertyname)
        classifer_struct.(propertyname{ii}) = classifer_obj.(propertyname{ii});
    end
    save('./additional_analyses/data/classifer_struct.mat', 'classifer_struct');
    
    propertyname = properties(noise_corr_obj);
    noise_corr_struct = struct(); %use a struct to store the information
    for ii = 1:numel(propertyname)
        noise_corr_struct.(propertyname{ii}) = noise_corr_obj.(propertyname{ii});
    end
    save('./additional_analyses/data/noise_corr_struct.mat', 'noise_corr_struct');
    
else
    load('./additional_analyses/data/fixation_struct.mat', 'fixation_struct');
    load('./additional_analyses/data/classifer_struct.mat', 'classifer_struct');
    load('./additional_analyses/data/noise_corr_struct.mat', 'noise_corr_struct');
    
    fixation_obj = fixation_session(fixation_struct);
    classifer_obj = novelty_classifier_class(classifer_struct);
    noise_corr_obj = noise_correlation_class('noveltyresponsive',noise_corr_struct);
    
end

fixation_obj.fixation_analysis(plotpath)

classifer_obj.novelty_classifier_analysis(plotpath);

noise_corr_obj.noise_correlation_analysis(plotpath, 'normalize', 1);

noise_corr_obj.noise_variance_analysis(plotpath);

toc






