% test the total number of neuronlist in the raw file and simplified files
% are matched

%Slayer
filepath = 'Y:\PLEXON_GRAYARRAY_Slayer\NovelFractalLearning_Slayer_for_share';
DDD = dir(fullfile(filepath,'*tasksession_mega_file.mat'));
%Lemmy
filepath = 'Y:\PLEXON_GRAYARRAY_LEMMY\NovelfractalLearning_Lemmy_for_share';
DDD = [DDD;dir(fullfile(filepath,'*tasksession_mega_file.mat'))];

Neuronlist_all_simp = [];
for i = 1:min(length(DDD),inf)
    filename = DDD(i).name;
    filepath = DDD(i).folder;
    
    load(fullfile(filepath, filename),'Neuronlist');
    Neuronlist_all_simp = [Neuronlist_all_simp, Neuronlist];
    
    fprintf('%d out of %d files done\n', i, numel(DDD));
end