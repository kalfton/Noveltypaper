%% 5 day's plot: Only Monkty S has data for it.
appearnum = 5;

variablenames = {'Learning_Familiar', 'Learning_Novel', 'Learning_1day', 'Learning_2day', 'Learning_3day', 'Learning_4day', 'Learning_5day'};
%Monkey S
fractalIDset = {6300:6303, 7999, [7410:7411,7420:7421], [7412,7422], [7413,7423], [7414,7424], [7415,7425]};


% choose the neurons in the session which has multiday fractal
logical_multiday = cellfun(@(x) ~isempty(x.('FR7410'))|~isempty(x.('FR7411')), {Neuronlist_all(:).learning})';

logical_multiday = cellfun(@(x) ~isempty(x.('FR7410')) | ~isempty(x.('FR7411')), {Neuronlist_all(:).learning})';
logical_multiday = logical_multiday & cellfun(@(x) (numel(x.('learningdate'))==5), {Neuronlist_all(:).learning})';

logical_for_neuronlist = {[Neuronlist_all(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_all(:).pred_nov_vs_fam]'>0 & logical_multiday};

Select_criteria = {'Novelty excited'};

Learning_5_day_struct = struct([]);
for xxx = 1: length(Neuronlist_all)
    familiar_FRstart = [];
    familiar_FRend = [];
    for i = 1:length(variablenames)
        % frac in learning trial
        fracID = fractalIDset{i};
        FRstart = []; %length = appearnum*length(fracID)
        FRend = [];
        for ij = 1:length(fracID)
            structname = ['FR', num2str(fracID(ij))];
            if isfield(Neuronlist_all(xxx).learning, structname)
                tempFR = Neuronlist_all(xxx).learning.(structname);
                if length(tempFR)>2*appearnum
                    FRstart = [FRstart; Neuronlist_all(xxx).learning.(structname)(1:appearnum)];
                    FRend = [FRend; Neuronlist_all(xxx).learning.(structname)(end-appearnum+1:end)];
                end
            end
        end
        %end
        
        Learning_5_day_struct(xxx).([variablenames{i} '_start']) = nanmean(FRstart);
        Learning_5_day_struct(xxx).([variablenames{i} '_end']) = nanmean(FRend);
        
        
        
        % ROC: novel vs familiar, learning vs familiar
        if contains(variablenames{i},'Familiar')
            familiar_FRstart = FRstart;
            familiar_FRend = FRend;
            Learning_5_day_struct(xxx).([variablenames{i} '_startroc']) = nan;
            Learning_5_day_struct(xxx).([variablenames{i} '_endroc']) = nan;
            Learning_5_day_struct(xxx).([variablenames{i} '_startp']) = nan;
            Learning_5_day_struct(xxx).([variablenames{i} '_endp']) = nan;
        else
            if ~isempty(familiar_FRstart) && ~isempty(FRstart)
                roc_start = rocarea3(familiar_FRstart, FRstart);
                roc_end = rocarea3(familiar_FRend, FRend);
                roc_start_p = ranksum(familiar_FRstart, FRstart);
                roc_end_p = ranksum(familiar_FRend, FRend);
                
                Learning_5_day_struct(xxx).([variablenames{i} '_startroc']) = roc_start;
                Learning_5_day_struct(xxx).([variablenames{i} '_endroc']) = roc_end;
                Learning_5_day_struct(xxx).([variablenames{i} '_startp']) = roc_start_p;
                Learning_5_day_struct(xxx).([variablenames{i} '_endp']) = roc_end_p;
                
            else
                Learning_5_day_struct(xxx).([variablenames{i} '_startroc']) = nan;
                Learning_5_day_struct(xxx).([variablenames{i} '_endroc']) = nan;
                Learning_5_day_struct(xxx).([variablenames{i} '_startp']) = nan;
                Learning_5_day_struct(xxx).([variablenames{i} '_endp']) = nan;
            end
        end
    end
end
    
for xy = 1: length(logical_for_neuronlist)
Neuronlist_learning = Neuronlist_all(logical_for_neuronlist{xy});
Neuronum = sum(logical_for_neuronlist{xy});

figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
plot_nrow = 2;
plot_ncol= 3;
plotplacesetx = {1,2};
plotplacesety = {1,1};


% summarize the Roc of five days:
nsubplot(plot_nrow, plot_ncol, plotplacesety{1}, plotplacesetx{1});
ROC_start_mean = [];
ROC_start_std = [];
ROC_end_mean = [];
ROC_end_std = [];
for i = 3: length(variablenames)
    All_ROC_start = [Learning_5_day_struct(logical_for_neuronlist{xy}).([variablenames{i} '_startroc'])];
    All_ROC_end = [Learning_5_day_struct(logical_for_neuronlist{xy}).([variablenames{i} '_endroc'])];
    
    ROC_start_mean = [ROC_start_mean, nanmean(All_ROC_start)];
    ROC_start_std = [ROC_start_std,nanstd(All_ROC_start)/sqrt(length(All_ROC_start))];
    ROC_end_mean = [ROC_end_mean, nanmean(All_ROC_end)];
    ROC_end_std = [ROC_end_std,nanstd(All_ROC_end)/sqrt(length(All_ROC_end))];
end


errorbar([1:5], ROC_start_mean, ROC_start_std, 'm');hold on;
errorbar([1:5]+0.5, ROC_end_mean, ROC_end_std, 'b');

% make the plots looks better

plot([[1:5];[1:5]+0.5], [ROC_start_mean; ROC_end_mean], 'k');
plot([[1:4]+0.5;[2:5]], [ROC_end_mean(1:4); ROC_start_mean(2:5)], 'color', [0.5,0.5,0]);

plot(get(gca, 'xlim'), [1,1]*0.5,'color', [0.8,0.8,0.8]);

xlim([0.5,7.5]);
ylim([0.5,0.65]);
legend({'Begin', 'End'}, 'Location','northeast');
legend('boxoff');
ylabel('ROC');
xlabel('Days')
set(gca, 'xtick', [1:5]);

print(gcf,'-dpdf', '-painters',[plotpath '/Learning_across_day_' Select_criteria{xy} '_5days_separate.pdf']);

end



