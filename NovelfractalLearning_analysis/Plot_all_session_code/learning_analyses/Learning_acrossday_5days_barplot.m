%% 5 day's plot: Only Slayer has data for it.
appearnum = 5;

variablenames = {'Learning_Familiar', 'Learning_Novel', 'Learning_1day', 'Learning_2day', 'Learning_3day', 'Learning_4day', 'Learning_5day'};
%Slayer
fractalIDset = {6300:6303, 7999, [7410:7411,7420:7421], [7412,7422], [7413,7423], [7414,7424], [7415,7425]};
fractaldateset = {nan, nan, nan, nan, nan, nan, nan};


% choose the neurons in the session which has multiday fractal
logical_multiday = cellfun(@(x) ~isempty(x.('FR7410'))|~isempty(x.('FR7411')), {Neuronlist_good(:).learning})';

logical_multiday = cellfun(@(x) ~isempty(x.('FR7410')) | ~isempty(x.('FR7411')), {Neuronlist_good(:).learning})';
logical_multiday = logical_multiday & cellfun(@(x) (numel(x.('learningdate'))==5), {Neuronlist_good(:).learning})';% || numel(x.('learningdate'))==1 && x.('learningdate')>1), {Neuronlist_good(:).learning})';
% Only include Slayer's neuron

logical_for_neuronlist = {[Neuronlist_good(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]'>0 & logical_multiday};

Select_criteria = {'Novelty excited'};


for xxx = 1: length(Neuronlist_good)
    familiar_FRstart = [];
    familiar_FRend = [];
    for i = 1:length(variablenames)
        % frac in learning trial
        fracID = fractalIDset{i};
        FRstart = []; %length = appearnum*length(fracID)
        FRend = [];
        %if any(isnan(fractaldateset{i})) || ismember(Neuronlist_good(xxx).learning.learningdate(1),fractaldateset{i})
        for ij = 1:length(fracID)
            structname = ['FR', num2str(fracID(ij))];
            if isfield(Neuronlist_good(xxx).learning, structname)
                tempFR = Neuronlist_good(xxx).learning.(structname);
                if length(tempFR)>2*appearnum
                    FRstart = [FRstart; Neuronlist_good(xxx).learning.(structname)(1:appearnum)];
                    FRend = [FRend; Neuronlist_good(xxx).learning.(structname)(end-appearnum+1:end)];
                end
            end
        end
        %end
        
        Neuronlist_good(xxx).([variablenames{i} '_start']) = nanmean(FRstart);
        Neuronlist_good(xxx).([variablenames{i} '_end']) = nanmean(FRend);
        
        
        
        % ROC: novel vs familiar, learning vs familiar
        if contains(variablenames{i},'Familiar')
            familiar_FRstart = FRstart;
            familiar_FRend = FRend;
            Neuronlist_good(xxx).([variablenames{i} '_startroc']) = nan;
            Neuronlist_good(xxx).([variablenames{i} '_endroc']) = nan;
            Neuronlist_good(xxx).([variablenames{i} '_startp']) = nan;
            Neuronlist_good(xxx).([variablenames{i} '_endp']) = nan;
        else
            if ~isempty(familiar_FRstart) && ~isempty(FRstart)
                roc_start = rocarea3(familiar_FRstart, FRstart);
                roc_end = rocarea3(familiar_FRend, FRend);
                roc_start_p = ranksum(familiar_FRstart, FRstart);
                roc_end_p = ranksum(familiar_FRend, FRend);
                
                Neuronlist_good(xxx).([variablenames{i} '_startroc']) = roc_start;
                Neuronlist_good(xxx).([variablenames{i} '_endroc']) = roc_end;
                Neuronlist_good(xxx).([variablenames{i} '_startp']) = roc_start_p;
                Neuronlist_good(xxx).([variablenames{i} '_endp']) = roc_end_p;
                
            else
                Neuronlist_good(xxx).([variablenames{i} '_startroc']) = nan;
                Neuronlist_good(xxx).([variablenames{i} '_endroc']) = nan;
                Neuronlist_good(xxx).([variablenames{i} '_startp']) = nan;
                Neuronlist_good(xxx).([variablenames{i} '_endp']) = nan;
            end
        end
        
        if contains(variablenames{i},'Novel')
            Novel_FRstart = FRstart;
            Novel_FRend = FRend;
            
        elseif contains(variablenames{i},'1day')
            Learning1day_FRstart = FRstart;
            Learning1day_FRend = FRend;
        elseif contains(variablenames{i},'2day')
            Learning2day_FRstart = FRstart;
            Learning2day_FRend = FRend;
        end
    end
    
    
    %%
    
    
    n = floor(numel(Learning1day_FRstart)/2);
    try
        Neuronlist_good(xxx).withindaylearningroc = rocarea3(Learning1day_FRend(1:n), Learning1day_FRstart(1:n));
        Neuronlist_good(xxx).withindaylearningroc_p = ranksum(Learning1day_FRend(1:n), Learning1day_FRstart(1:n));
        
    catch
        Neuronlist_good(xxx).withindaylearningroc = nan;
        Neuronlist_good(xxx).withindaylearningroc_p = nan;
    end
    
    try
        Neuronlist_good(xxx).acrossdayforgetroc = rocarea3(Learning1day_FRend((n+1):end), Learning2day_FRstart);
        Neuronlist_good(xxx).acrossdayforgetroc_p = ranksum(Learning1day_FRend(1:n), Learning2day_FRstart(1:n));
        
    catch
        Neuronlist_good(xxx).acrossdayforgetroc = nan;
        Neuronlist_good(xxx).acrossdayforgetroc_p = nan;
    end
    
    %%
    Neuronlist_good(xxx).learning_ind = max(min((Neuronlist_good(xxx).Learning_Novel_end - Neuronlist_good(xxx).Learning_1day_end)...
        /(Neuronlist_good(xxx).Learning_Novel_end - Neuronlist_good(xxx).Learning_Familiar_end),1),0);
    
end

novelty_resp_logic = [Neuronlist_good(:).withindaylearningroc_p]<StatisticalThreshold & [Neuronlist_good(:).withindaylearningroc]>0;

    
for xy = 1: length(logical_for_neuronlist)
Neuronlist_learning = Neuronlist_good(logical_for_neuronlist{xy});
Neuronum = sum(logical_for_neuronlist{xy});

for i = 1:length(variablenames)
    eval([variablenames{i} '_start = [Neuronlist_learning(:).' variablenames{i} '_start]'';']);
    eval([variablenames{i} '_end = [Neuronlist_learning(:).' variablenames{i} '_end]'';']);
    eval([variablenames{i} '_startroc = [Neuronlist_learning(:).' variablenames{i} '_startroc]'';']);
    eval([variablenames{i} '_endroc = [Neuronlist_learning(:).' variablenames{i} '_endroc]'';']);
    eval([variablenames{i} '_startp = [Neuronlist_learning(:).' variablenames{i} '_startp]'';']);
    eval([variablenames{i} '_endp = [Neuronlist_learning(:).' variablenames{i} '_endp]'';']);
    
end

learning_ind = [Neuronlist_learning(:).learning_ind]';




figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
plotplacesetx = {11:45,11:45};
plotplacesety = {121:169, 1:49};
% Write down the criteria
nsubplot(300,210, plotplacesety{end}, plotplacesetx{end});
text(1,1,Select_criteria{xy}, 'Fontsize', 10);
text(1,0,['N = ' num2str(Neuronum)], 'Fontsize', 10);
text(1,0.2,['Frac N = ' num2str(appearnum)], 'Fontsize', 10);

axis off;


% summarize the Roc of five days:
nsubplot(300,210, plotplacesety{1}, plotplacesetx{1});
ROC_start_mean = [];
ROC_start_std = [];
ROC_end_mean = [];
ROC_end_std = [];
for i = 3: length(variablenames)
    eval(['All_ROC_start =' variablenames{i} '_startroc;']);
    eval(['All_ROC_end =' variablenames{i} '_endroc;']);
    
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

plot(get(gca, 'xlim'), [1,1]*nanmean(Learning_Novel_startroc),'m');
plot(get(gca, 'xlim'), [1,1]*nanmean(Learning_Novel_endroc),'b');
plot(get(gca, 'xlim'), [1,1]*0.5,'color', [0.8,0.8,0.8]);

xlim([0.5,7.5]);
ylim([0.5,0.65]);
legend({'Begin', 'End', 'withinday', 'acrossday'}, 'Location','northeast');
legend('boxoff');
ylabel('ROC');
xlabel('Days')
set(gca, 'xtick', [1:5]);

print(gcf,'-dpdf', '-painters',[plotpath '/Learning_across_day_' Select_criteria{xy} '_5days_separate.pdf']);

end



