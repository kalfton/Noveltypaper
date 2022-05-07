%% Across day bar plots
% choose the neurons in the session which has multiday fractal
logical_multiday = cellfun(@(x) ~isempty(x.('FR7410')) | ~isempty(x.('FR7411')), {Neuronlist_all(:).learning})';
logical_multiday = logical_multiday & cellfun(@(x) (numel(x.('learningdate'))==5 || numel(x.('learningdate'))==1 && x.('learningdate')>1), {Neuronlist_all(:).learning})';
logical_multiday = logical_multiday;

logical_for_neuronlist = {[Neuronlist_all(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_all(:).pred_nov_vs_fam]'>0 & logical_multiday
    };
    %[Neuronlist_all(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_all(:).pred_nov_vs_fam]'<0 & logical_multiday
Select_criteria = {'Novelty excited', 'Novelty inhibited'};
variablenames = {'Learning_1day', 'Learning_2day'};
appearnum = 5;

for xy = 1: length(logical_for_neuronlist)
Neuronlist_learning = Neuronlist_all(logical_for_neuronlist{xy});
Neuronum = sum(logical_for_neuronlist{xy});

for i = 1:length(variablenames)
    eval([variablenames{i} '_start = [Neuronlist_learning(:).' variablenames{i} '_start]'';']);
    eval([variablenames{i} '_end = [Neuronlist_learning(:).' variablenames{i} '_end]'';']);
    eval([variablenames{i} '_startroc = [Neuronlist_learning(:).' variablenames{i} '_startroc]'';']);
    eval([variablenames{i} '_endroc = [Neuronlist_learning(:).' variablenames{i} '_endroc]'';']);
    eval([variablenames{i} '_startp = [Neuronlist_learning(:).' variablenames{i} '_startp]'';']);
    eval([variablenames{i} '_endp = [Neuronlist_learning(:).' variablenames{i} '_endp]'';']);
    
end


figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
plot_nrow = 2;
plot_ncol= 3;
plotplacesety = {1,1,1};
plotplacesetx = {1,2,3};
% Write down the criteria etc.
nsubplot(plot_nrow, plot_ncol, plotplacesety{end}, plotplacesetx{end});
text(1,1,Select_criteria{xy}, 'Fontsize', 10);
text(1,0,['N = ' num2str(Neuronum)], 'Fontsize', 10);
text(1,0.2,['Frac N = ' num2str(appearnum)], 'Fontsize', 10);

axis off;


for i = 1:length(variablenames)
    % ROC
    if ~contains(variablenames{i},{'Familiar', 'Novel'}) %only plot learning
        eval(['All_ROC_start =' variablenames{i} '_startroc;']);
        eval(['All_ROC_end =' variablenames{i} '_endroc;']);
        nsubplot(plot_nrow, plot_ncol, plotplacesety{i}, plotplacesetx{i});
        plot([0.5,2.5], [0.5,0.5], 'color', [0.8,0.8,0.8])
        
        Datamean = [nanmean(All_ROC_start), nanmean(All_ROC_end)];
        Datastd = [nanstd(All_ROC_start)/sqrt(length(All_ROC_start)), nanstd(All_ROC_end)/sqrt(length(All_ROC_end))];
        [p_paired, ~] = signrank(All_ROC_start-All_ROC_end);
        title([variablenames{i} 'Fractal, p = ' mat2str(p_paired,3)], 'interpreter', 'none');
        
        % p ranksum between end of Ln day and begin of Ln+1 day
        if contains(variablenames{i},'Learning') && contains(variablenames{i},'day') && ~contains(variablenames{i},'1day')
            eval(['Prev_day_ROC_end = ' variablenames{i-1} '_endroc;']);
            try
                p_paired = signrank(All_ROC_start-Prev_day_ROC_end);
            catch
                p_paired = nan;
            end
            text(1,0.5, ['p preday = ' mat2str(p_paired,3)]);
            
            eval(['Prev_day_ROC_start = ' variablenames{i-1} '_startroc;']);
            try
                p = signrank(All_ROC_start-Prev_day_ROC_start);
            catch
                p = nan;
            end
            text(1,0.7, ['p day1vsday2 start=' mat2str(p,3)]);
        end
        
        p_signrank = signrank(All_ROC_end,0.5);
        text(2,0.5, ['p signrank = ' mat2str(p_signrank,3)]);
        
        errorbar([1,2],Datamean, Datastd, 'r');
        set(gca, 'xtick', [1,2], 'xticklabel', {'Begin', 'End'});
        xlim([0.5,2.5]);
        ylim([0.4,0.7]);
        ylabel('ROC');
    end
end


print(gcf,'-dpdf', '-painters',[plotpath '/Learning_across_day_' Select_criteria{xy} '.pdf']);

end


%% surprise effect in learning
figure;
novelty_resp_logic = [Neuronlist_all(:).P_pred_nov_vs_fam]<StatisticalThreshold & [Neuronlist_all(:).pred_nov_vs_fam]>0 ...
    & ~isnan([Neuronlist_all(:).learning_surprise])...
    & logical_multiday';
title(['Learning trial type testing surprise, n = ' mat2str(sum(novelty_resp_logic))]);
learning_surprise = [Neuronlist_all(novelty_resp_logic).learning_surprise];


Datamean = [mean(learning_surprise)];
Datastd = [std(learning_surprise)/sqrt(length(learning_surprise))];

bar([1],Datamean, 'b'); hold on;
errorbar([1],Datamean, Datastd, 'r.');

set(gca, 'xtick', [1], 'xticklabel', {'Learning surprise'});
ylabel('indices');
xlim([0,2]);

y_limit = get(gca, 'ylim');
p = signrank(learning_surprise);
text(1,y_limit(1)*0.7,['p = ' mat2str(p,3)]);
print(gcf,'-dpdf', '-painters',[plotpath '/Learning_surprise.pdf']);