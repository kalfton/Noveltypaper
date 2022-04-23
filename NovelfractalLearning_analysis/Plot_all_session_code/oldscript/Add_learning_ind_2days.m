function Neuronlist_good = Add_learning_ind_2days(Neuronlist_good, pars)
% this script add learning index to neuronlist

%close all;
%shuffling_num = 10000;
run_longtimebootstrapping = 0;
StatisticalThreshold = 0.01;
% 
%Slayer
fractalIDset_Slayer = {6300:6303, 7999, [7410:7411,7420:7421], [7412,7422,7413,7423, 7414,7424, 7415,7425]};
%logical_Slayer = cellfun(@(x) (strcmpi(x ,'Slayer')), {Neuronlist_good(:).monkey})';
%Lemmy
fractalIDset_Lemmy = {6300:6307, 7999, [7300:7307], [7410:7411,7420:7421]};
%logical_Lemmy = cellfun(@(x) (strcmpi(x ,'Lemmy')), {Neuronlist_good(:).monkey})';

% plot the neuron's learning curve separately for novelty excited,
% inhibited, and other neuron.
% And recency
% choose the neurons in the sessions which has multiday fractal
logical_multiday = cellfun(@(x) ~isempty(x.('FR7410')) | ~isempty(x.('FR7411')), {Neuronlist_good(:).learning})';
logical_multiday = logical_multiday & cellfun(@(x) (numel(x.('learningdate'))==5 || numel(x.('learningdate'))==1 && x.('learningdate')>1), {Neuronlist_good(:).learning})';
%logical_multiday = logical_multiday;


appearnum = 5;

variablenames = {'Learning_Familiar', 'Learning_Novel', 'Learning_1day', 'Learning_2day'};
fractaldateset = {nan, nan, nan, nan, nan, nan, nan};



for xxx = 1: length(Neuronlist_good)
    %%
    
    familiar_FRstart = [];
    familiar_FRend = [];
    
    if strcmpi(Neuronlist_good(xxx).monkey, 'Lemmy')
        fractalIDset = fractalIDset_Lemmy;
    elseif strcmpi(Neuronlist_good(xxx).monkey, 'Slayer')
        fractalIDset = fractalIDset_Slayer;
    end
    for i = 1:length(variablenames)
        % frac in learning trial
        fracID = fractalIDset{i};
        FRstart = []; %length = appearnum*length(fracID)
        FRend = [];
        if any(isnan(fractaldateset{i})) || ismember(Neuronlist_good(xxx).learning.learningdate(1),fractaldateset{i})
            for ij = 1:length(fracID)
                structname = ['FR', num2str(fracID(ij))];
                if isfield(Neuronlist_good(xxx).learning, structname)
                    tempFR = Neuronlist_good(xxx).learning.(structname);
                    if length(tempFR)>=2*appearnum
                        FRstart = [FRstart; Neuronlist_good(xxx).learning.(structname)(1:appearnum)];
                        FRend = [FRend; Neuronlist_good(xxx).learning.(structname)(end-appearnum+1:end)];
                    end
                end
            end
        end
        
        Neuronlist_good(xxx).([variablenames{i} '_start']) = nanmean(FRstart);
        Neuronlist_good(xxx).([variablenames{i} '_end']) = nanmean(FRend);
        
        Neuronlist_good(xxx).([variablenames{i} '_start_all']) = FRstart;
        Neuronlist_good(xxx).([variablenames{i} '_end_all']) = FRend;
        
        % ROC: novel vs familiar, learning vs familiar
        % initiate
        roc_start = nan;
        roc_end = nan;
        roc_start_p = nan;
        roc_end_p = nan;
        
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
            else
                roc_start = nan;
                roc_end = nan;
                roc_start_p = nan;
                roc_end_p = nan;
            end
            Neuronlist_good(xxx).([variablenames{i} '_startroc']) = roc_start;
            Neuronlist_good(xxx).([variablenames{i} '_endroc']) = roc_end;
            Neuronlist_good(xxx).([variablenames{i} '_startp']) = roc_start_p;
            Neuronlist_good(xxx).([variablenames{i} '_endp']) = roc_end_p;
        end
        if contains(variablenames{i},'Novel')
            Novel_FRstart = FRstart;
            Novel_FRend = FRend;
            
        elseif contains(variablenames{i},'1day')
            Learning1day_FRstart = FRstart;
            Learning1day_FRend = FRend;
            Learning1day_startroc = roc_start;
            Learning1day_endroc = roc_end;
        elseif contains(variablenames{i},'2day')
            Learning2day_FRstart = FRstart;
            Learning2day_FRend = FRend;
            Learning2day_startroc = roc_start;
            Learning2day_endroc = roc_end;
        end
    end
    start_FR_NonLearning_Novel = [];
    end_FR_NonLearning_Novel = [];
    tempFR = Neuronlist_good(xxx).nonlearning.FR7999;
    % old converter's nonlearning are not z-scored
%     if Neuronlist_good(xxx).zscorestd>1e-7
%         tempFR = (tempFR-Neuronlist_good(xxx).zscoremean)/Neuronlist_good(xxx).zscorestd;
%     else
%         tempFR = (tempFR-Neuronlist_good(xxx).zscoremean);
%     end
    if length(tempFR)>2*appearnum
        start_FR_NonLearning_Novel = tempFR(1:appearnum);
        end_FR_NonLearning_Novel = tempFR(end-appearnum+1:end);
    end
%     start_FR_Learning_Novel = Neuronlist_good(xxx).Learning_Novel_start;
%     end_FR_Learning_Novel = Neuronlist_good(xxx).Learning_Novel_end;
%     start_FR_Learning_Familiar = Neuronlist_good(xxx).Learning_Familiar_start;
%     end_FR_Learning_Familiar = Neuronlist_good(xxx).Learning_Familiar_end;
    Neuronlist_good(xxx).nonLearning_Novel_start = mean(start_FR_NonLearning_Novel);
    Neuronlist_good(xxx).nonLearning_Novel_end = mean(end_FR_NonLearning_Novel);
    
    try
        Neuronlist_good(xxx).learning_surprise = rocarea3(Novel_FRstart, Novel_FRend) - rocarea3(start_FR_NonLearning_Novel, end_FR_NonLearning_Novel);
        Neuronlist_good(xxx).learning_recency  = rocarea3(Learning1day_FRend, Learning2day_FRstart) - rocarea3(end_FR_NonLearning_Novel, start_FR_NonLearning_Novel);
        Neuronlist_good(xxx).learning_recency_LvsF = Learning2day_startroc-Learning1day_endroc;
    catch
        Neuronlist_good(xxx).learning_surprise = nan;
        Neuronlist_good(xxx).learning_recency  = nan;
        Neuronlist_good(xxx).learning_recency_LvsF = nan;
    end
    
    %rocarea3(Learning1day_FRend, Learning2day_FRstart) - rocarea3(end_FR_NonLearning_Novel, start_FR_NonLearning_Novel);
    %rocarea3(familiar_FRend, familiar_FRstart) - rocarea3(end_FR_NonLearning_Novel, start_FR_NonLearning_Novel);
    if isempty(Neuronlist_good(xxx).learning_surprise)
        Neuronlist_good(xxx).learning_surprise = nan;
    end
    if isempty(Neuronlist_good(xxx).learning_recency)
        Neuronlist_good(xxx).learning_recency = nan;
    end
    %
    Learning_Novel_all_FR = Neuronlist_good(xxx).learning.FR7999;
    Learning_Familiar_all_FR = [Neuronlist_good(xxx).learning.FR6300
        Neuronlist_good(xxx).learning.FR6301
        Neuronlist_good(xxx).learning.FR6302
        Neuronlist_good(xxx).learning.FR6303];
    Neuronlist_good(xxx).learning_nov_fam_roc = rocarea3(Learning_Familiar_all_FR, Learning_Novel_all_FR);
    Neuronlist_good(xxx).learning_nov_fam_p = ranksum(Learning_Familiar_all_FR, Learning_Novel_all_FR);
    %
    
    n = floor(numel(Learning1day_FRstart)/2);
    n2 = floor(numel(familiar_FRstart)/5/2)*5;
    try
        Neuronlist_good(xxx).withindaylearningroc = rocarea3(Learning1day_FRend(1:n), Learning1day_FRstart(1:n))-rocarea3(familiar_FRend(1:n2), familiar_FRstart(1:n2))+0.5;
        Neuronlist_good(xxx).withindaylearningroc_p = ranksum(Learning1day_FRend(1:n), Learning1day_FRstart(1:n));
        if isempty(Neuronlist_good(xxx).withindaylearningroc)
            Neuronlist_good(xxx).withindaylearningroc = nan;
        end
        
        %         Neuronlist_good(xxx).withindaylearningroc = (rocarea3(Learning1day_FRend(1:n), Learning1day_FRstart(1:n))-0.5)/(rocarea3(Learning_Familiar_all_FR, Learning_Novel_all_FR)-0.5)+0.5;
        %         Neuronlist_good(xxx).withindaylearningroc_p = ranksum(Learning1day_FRend(1:n), Learning1day_FRstart(1:n));
        %
        %         if Neuronlist_good(xxx).withindaylearningroc>1
        %             Neuronlist_good(xxx).withindaylearningroc = 1;
        %         elseif Neuronlist_good(xxx).withindaylearningroc<0
        %             Neuronlist_good(xxx).withindaylearningroc = 0;
        %         end
        %
    catch
        Neuronlist_good(xxx).withindaylearningroc = nan;
        Neuronlist_good(xxx).withindaylearningroc_p = nan;
    end
    
    try
        Neuronlist_good(xxx).acrossdayforgetroc = rocarea3(Learning1day_FRend((n+1):end), Learning2day_FRstart)-rocarea3(familiar_FRend((n2+1):end), familiar_FRstart((n2+1):end))+0.5;
        Neuronlist_good(xxx).acrossdayforgetroc_p = ranksum(Learning1day_FRend((n+1):end), Learning2day_FRstart);
        if isempty(Neuronlist_good(xxx).acrossdayforgetroc)
            Neuronlist_good(xxx).acrossdayforgetroc = nan;
        end
        
        %         Neuronlist_good(xxx).acrossdayforgetroc = (rocarea3(Learning1day_FRend((n+1):end), Learning2day_FRstart)-0.5)/(rocarea3(Learning1day_FRend(1:n), Learning1day_FRstart(1:n))-0.5)+0.5;
        %         Neuronlist_good(xxx).acrossdayforgetroc_p = ranksum(Learning1day_FRend(1:n), Learning2day_FRstart(1:n));
        %
        %         if Neuronlist_good(xxx).acrossdayforgetroc>1
        %             Neuronlist_good(xxx).acrossdayforgetroc = 1;
        %         elseif Neuronlist_good(xxx).acrossdayforgetroc<0
        %             Neuronlist_good(xxx).acrossdayforgetroc = 0;
        %         end
        
    catch
        Neuronlist_good(xxx).acrossdayforgetroc = nan;
        Neuronlist_good(xxx).acrossdayforgetroc_p = nan;
    end
    
    %% New way to do learning withinday and forgetting across day index
    try
        %averaged variance of neuron's firing rate
        sigma2 = var(Novel_FRstart,1)*numel(Novel_FRstart)+var(Novel_FRend,1)*numel(Novel_FRend)+...
            var(familiar_FRend,1)*numel(familiar_FRend)+var(familiar_FRstart,1)*numel(familiar_FRstart);
        sigma2 = sigma2/(numel(Novel_FRstart)+numel(Novel_FRend)+numel(familiar_FRend)+numel(familiar_FRstart)-4);
        
        mean_Novel_start = mean(Novel_FRstart);
        mean_Familiar_start = mean(familiar_FRstart);
        mean_Novel_end = mean(Novel_FRend);
        mean_Familiar_end = mean(familiar_FRend);
        
        % classifying learning fractal's firing rate between 0,1(Familiar/Novel, with assumption of Gaussian distribution. Actually is
        % using a sigmoid function
        
        sigmoidfunstart = @(x) 1./(1+exp(-(mean_Novel_start-mean_Familiar_start)./sigma2*(x-(mean_Novel_start+mean_Familiar_start)/2)));
        sigmoidfunend = @(x) 1./(1+exp(-(mean_Novel_end-mean_Familiar_end)./sigma2*(x-(mean_Novel_end+mean_Familiar_end)/2)));
        
        n = floor(numel(Learning1day_FRend)/2);
        normal_Learning1day_endmean = mean(sigmoidfunend(Learning1day_FRend));
        normal_Learning1day_endmean1 = mean(sigmoidfunend(Learning1day_FRend(1:n)));
        normal_Learning1day_endmean2 = mean(sigmoidfunend(Learning1day_FRend((n+1):end)));
        normal_Learning1day_startmean = mean(sigmoidfunstart(Learning1day_FRstart));
        normal_Learning2day_startmean = mean(sigmoidfunstart(Learning2day_FRstart));
        normal_Novel_startmean = mean(sigmoidfunstart(Novel_FRstart));
        
        Neuronlist_good(xxx).withindaylearning_newindex = normal_Learning1day_startmean - normal_Learning1day_endmean1;
        Neuronlist_good(xxx).acrossdayforget_newindex = (normal_Learning2day_startmean - normal_Learning1day_endmean2);%/Neuronlist_good(xxx).withindaylearning_newindex;
        
        Neuronlist_good(xxx).withindaylearning_newindex_p = ranksum(sigmoidfunstart(Learning1day_FRstart), sigmoidfunend(Learning1day_FRend(1:n)));
        Neuronlist_good(xxx).acrossdayforget_newindex_p = ranksum(sigmoidfunstart(Learning2day_FRstart), sigmoidfunend(Learning1day_FRend((n+1):end)));
        
        Neuronlist_good(xxx).normal_Learning_1day_start = mean(sigmoidfunstart(Learning1day_FRstart));
        Neuronlist_good(xxx).normal_Learning_1day_end = mean(sigmoidfunend(Learning1day_FRend));
        Neuronlist_good(xxx).normal_Learning_2day_start = mean(sigmoidfunstart(Learning2day_FRstart));
        Neuronlist_good(xxx).normal_Learning_2day_end = mean(sigmoidfunend(Learning2day_FRend));
        Neuronlist_good(xxx).normal_Learning_Novel_start = mean(sigmoidfunstart(Novel_FRstart));
        Neuronlist_good(xxx).normal_Learning_Novel_end = mean(sigmoidfunend(Novel_FRend));
        Neuronlist_good(xxx).normal_Learning_Familiar_start = mean(sigmoidfunstart(familiar_FRstart));
        Neuronlist_good(xxx).normal_Learning_Familiar_end = mean(sigmoidfunend(familiar_FRend));
        
    catch
        Neuronlist_good(xxx).withindaylearning_newindex = nan;
        Neuronlist_good(xxx).acrossdayforget_newindex = nan;
        Neuronlist_good(xxx).withindaylearning_newindex_p = nan;
        Neuronlist_good(xxx).acrossdayforget_newindex_p = nan;
        
        Neuronlist_good(xxx).normal_Learning_1day_start = nan;
        Neuronlist_good(xxx).normal_Learning_1day_end = nan;
        Neuronlist_good(xxx).normal_Learning_2day_start = nan;
        Neuronlist_good(xxx).normal_Learning_2day_end = nan;
        Neuronlist_good(xxx).normal_Learning_Novel_start = nan;
        Neuronlist_good(xxx).normal_Learning_Novel_end = nan;
        Neuronlist_good(xxx).normal_Learning_Familiar_start = nan;
        Neuronlist_good(xxx).normal_Learning_Familiar_end = nan;
    end
    
    %% other old indices
    start_roc_1day = Neuronlist_good(xxx).Learning_1day_startroc;
    end_roc_1day = Neuronlist_good(xxx).Learning_1day_endroc;
    
    start_roc_2day = Neuronlist_good(xxx).Learning_2day_startroc;
    end_roc_2day = Neuronlist_good(xxx).Learning_2day_endroc;
    
    try
        Neuronlist_good(xxx).Learning_recency_ind = (end_roc_2day - start_roc_2day)/(end_roc_1day - start_roc_1day);
    catch
        Neuronlist_good(xxx).Learning_recency_ind = nan;
    end
    
    learning_ind = (Neuronlist_good(xxx).Learning_Novel_end - Neuronlist_good(xxx).Learning_1day_end)...
        /(Neuronlist_good(xxx).Learning_Novel_end - Neuronlist_good(xxx).Learning_Familiar_end);
    if ~isnan(learning_ind)
        learning_ind = max(min(learning_ind,1),0);
    end
    Neuronlist_good(xxx).learning_ind = learning_ind;
    %%
end

if pars.subtractbysession
    % to avoid session by session variance. Subtract the mean of each
    % session
    novelty_resp_logic = [Neuronlist_good(:).P_pred_nov_vs_fam]<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]>0 & logical_multiday'...
        & ~isnan([Neuronlist_good(:).withindaylearning_newindex]) & ~isnan([Neuronlist_good(:).acrossdayforget_newindex]);
    
    if sum(novelty_resp_logic)>=5 % only include sessions that has more than 5 novelty selective neurons.
        mean_withindaylearning = mean([Neuronlist_good(novelty_resp_logic).withindaylearning_newindex]);
        mean_acrossdayforget = mean([Neuronlist_good(novelty_resp_logic).acrossdayforget_newindex]);
        for xxx = 1: length(Neuronlist_good)
            Neuronlist_good(xxx).withindaylearning_newindex = Neuronlist_good(xxx).withindaylearning_newindex - mean_withindaylearning;
            Neuronlist_good(xxx).acrossdayforget_newindex = Neuronlist_good(xxx).acrossdayforget_newindex - mean_acrossdayforget;
            if novelty_resp_logic(xxx)
                Neuronlist_good(xxx).learningforgetinganalysis = true;
            else
                Neuronlist_good(xxx).learningforgetinganalysis = false;
            end
        end
    else
        for xxx = 1: length(Neuronlist_good)
            Neuronlist_good(xxx).learningforgetinganalysis = false;
        end
    end
    
else
    novelty_resp_logic = [Neuronlist_good(:).P_pred_nov_vs_fam]<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]>0 & logical_multiday'...
        & ~isnan([Neuronlist_good(:).withindaylearning_newindex]) & ~isnan([Neuronlist_good(:).acrossdayforget_newindex]);
    
    for xxx = 1: length(Neuronlist_good)
        if novelty_resp_logic(xxx)
            Neuronlist_good(xxx).learningforgetinganalysis = true;
        else
            Neuronlist_good(xxx).learningforgetinganalysis = false;
        end
    end
    
end

end