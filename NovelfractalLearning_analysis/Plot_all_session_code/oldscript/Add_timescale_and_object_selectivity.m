function Neuronlist = Add_timescale_and_object_selectivity(Neuronlist, pars)
% pars.filename; pars.pathname


addpath('Y:\Kaining\NovelfractalLearning_analysis\code for parcluster\RIS_parallel\utils')
rng(0);
load(fullfile(pars.pathname, pars.filename),'Generaltask');

% Channelnames = [Generaltask.channelname];
% if numel(Channelnames)~=numel(Neuronlist)
%     error('Neuronlist and Channelnames mismatch');
% end


trialsplitlogical = true(size(Generaltask.successtrial));
fracsplitlogical = kron(trialsplitlogical',true(1,3));

successtrial_split = find(Generaltask.successtrial == 1);

%% generate fractal sets
frac_sets = make_frac_sets(Generaltask, fracsplitlogical);

%% loop through channels (only for SPK)
for xxx = 1:length(Neuronlist)
    channelname = Neuronlist(xxx).name;
    
    if contains(channelname, 'MUA')% skip caluculating MUA to save time.
        Neuronlist(xxx).intrinsic_time = nan;
        Neuronlist(xxx).intrinsic_A = nan;
        Neuronlist(xxx).intrinsic_B = nan;
        
        Neuronlist(xxx).P_object_select_novel_start = nan;
        Neuronlist(xxx).obj_select_ind_novel_start = nan;
        Neuronlist(xxx).P_object_select_novel_end = nan;
        Neuronlist(xxx).obj_select_ind_novel_end = nan;
        
        Neuronlist(xxx).P_object_select_recent = nan;
        Neuronlist(xxx).obj_select_ind_recent = nan;
        Neuronlist(xxx).P_object_select_nonrecent = nan;
        Neuronlist(xxx).obj_select_ind_nonrecent = nan;
        
        Neuronlist(xxx).P_object_select_surprise = nan;
        Neuronlist(xxx).obj_select_ind_surprise = nan;
        Neuronlist(xxx).P_object_select_nonsurprise = nan;
        Neuronlist(xxx).obj_select_ind_nonsurprise = nan;
        
        continue;
    else
        %load the channel
        load(fullfile(pars.pathname, pars.filename),channelname);
    end
    disp(channelname)
    
    
%     channeltype = [];
%     if contains(channelname,{'SPK','MUA'})
%         channeltype = 1; % 1 means spike channel
%         xxxy = find(strcmpi(channelname,{Neuronlist(:).name}));
%         Neuronlist(xxx) = Neuronlist(xxxy);
%     elseif contains(channelname,{'FP'})
%         channeltype = 2; % 2 means LFP channel
%         xxxy = find(strcmpi(channelname,{LFPlist(:).name}));
%         Neuronlist(xxx) = LFPlist(xxxy);
%     elseif contains(channelname,{'pupil'})
%         channeltype = 3; % 3 means pupil channel
%         xxxy = 1;
%         Neuronlist(xxx) = struct();
%         Neuronlist(xxx).name = 'pupil';
%     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Calculate neuron's PSTH and firing rate before fixation point
    binsize = 50;%ms
    
    Hist_count = [];
    %trial rasters
    eval(['sptimes = ' channelname '.sptimes;']);
    for x=1:size(Generaltask.trialstart,1)
        Hist_template = -700:binsize:0;
        Hist_raw = histcounts((sptimes{1,x}-Generaltask.trialstart(x,1)-Generaltask.timetargeton(x,1))*1000, Hist_template);
        Hist_count = [Hist_count; Hist_raw];
    end
    %Hist_sdf=plot_mean_psth({Hist_count},gauswindow_ms,1,size(Hist_count,2),1);
    Hist_template = Hist_template(1:end-1);
    
    
    
    %% Intrinsic timescale
    % autocorrelaton:
    autocorr_matrix = nan(numel(Hist_template)-1);
    for i= 1: numel(Hist_template)-1
        for j=1: numel(Hist_template)-i
            autocorr_matrix(i,j) = corr(Hist_count(successtrial_split,i), Hist_count(successtrial_split,i+j), 'Type','Pearson');
        end
    end
    
    % fit by exponiential function
    fitoptions = optimset('MaxIter',10000,'MaxFunEvals',20000);
    paras_init = [0.2,0,200];
    LB = [0, -Inf, 0];
    UB = [Inf, Inf, Inf];
    
    ave_cor = nanmean(autocorr_matrix,1);
    if ~any(isnan(ave_cor)) % check there are data point at each time step
        
        % to avoid fitting the refraction period compare the autocorrelation within 200ms, and use the highest one as start time interval
        [~,ind] = max(ave_cor(1:4));
        autocorr_matrix(:,1:ind-1) = nan;
        
        lossfun = @(x) Loss_exp_fun_cor(autocorr_matrix, x, binsize);
        exp_paras = fminsearchbnd(lossfun,paras_init,LB,UB,fitoptions);%%%
        
        Neuronlist(xxx).intrinsic_time = exp_paras(3);
        Neuronlist(xxx).intrinsic_A = exp_paras(1);
        Neuronlist(xxx).intrinsic_B = exp_paras(2);
        
        % Plot the fitted curve to check sanity
%         close all;
%         figure;hold on;
%         xdata = ones(size(autocorr_matrix,1),1)*binsize*(1:size(autocorr_matrix,1));
%         xdata = xdata(:);
%         ydata = autocorr_matrix(:);
%         scatter(xdata, ydata, 'k.');
%         xdata = binsize*(1:size(autocorr_matrix,1));
%         ydata = nanmean(autocorr_matrix,1);
%         plot(xdata,ydata, 'k');
%         exp_fun = @(k, A, B, t) A*(exp(-k*binsize/t)+B);
%         plot(binsize*(1:size(autocorr_matrix)), exp_fun(1:size(autocorr_matrix), exp_paras(1), exp_paras(2), exp_paras(3)));
%         % 
    else
        Neuronlist(xxx).intrinsic_time = nan;
        Neuronlist(xxx).intrinsic_A = nan;
        Neuronlist(xxx).intrinsic_B = nan;
    end
    
    
    
    All_Fractal_FR = Neuronlist(xxx).All_Fractal_FR;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % object selectivity & position effect: trial type 4 & 5
%     familiar45type = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6500:6510])));
%     FiringR_45 = All_Fractal_FR(familiar45type,1);
% 
%     [FR_subtracted, anova_result]= regression_anova(FiringR_45, Generaltask.Fractals(:,familiar45type), {'position','fractalID'});
%     
%     anova_p = anova_result.anova_p;
%     anova_table = anova_result.anova_table;
%     anova_stats = anova_result.anova_stats;
%     
%     % object selectivity index
%     var_residual = anova_table{4,5};
%     var_obj_ID = anova_table{3,5};
%     obj_select_ind = var_obj_ID/(var_residual+var_obj_ID); % object selectivity index is between 0 and 1, the higher value, the more selective it is.
%     Neuronlist(xxx).P_object_select_anova = anova_p(2);
%     Neuronlist(xxx).obj_select_ind = obj_select_ind;
    
    %% Novelty - object selectivity interaction
    
    % first one third of the repeating novelty vs. last one third of the
    % repeating novel fractals
    if strcmpi(Neuronlist(xxx).monkey, 'Slayer')
        learning1day_frac = [7410:7411,7420:7421];
    elseif strcmpi(Neuronlist(xxx).monkey, 'Lemmy')
        learning1day_frac = [7300:7307];
    else
        error('Monkey name invalid');
    end
    
    repeatingnovel_frac = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), learning1day_frac)));
    repeatingnovel_start = repeatingnovel_frac(1 : ceil(numel(repeatingnovel_frac)/3));
    repeatingnovel_end = repeatingnovel_frac(floor(numel(repeatingnovel_frac)*2/3) : end);
    FiringR_learning = All_Fractal_FR(repeatingnovel_frac,1);
    
    [~, anova_result]= regression_anova(All_Fractal_FR(repeatingnovel_start,1), Generaltask.Fractals(:,repeatingnovel_start), {'position','fractalID'});
    anova_p = anova_result.anova_p;
    anova_table = anova_result.anova_table;
    
    var_residual = anova_table{4,5};
    var_obj_ID = anova_table{3,5};
    obj_select_ind = var_obj_ID/(var_residual+var_obj_ID); % object selectivity index is between 0 and 1, the higher value, the more selective it is.
    Neuronlist(xxx).P_object_select_novel_start = anova_p(2);
    Neuronlist(xxx).obj_select_ind_novel_start = obj_select_ind;
    
    [~, anova_result]= regression_anova(All_Fractal_FR(repeatingnovel_end,1), Generaltask.Fractals(:,repeatingnovel_end), {'position','fractalID'});
    anova_p = anova_result.anova_p;
    anova_table = anova_result.anova_table;
    
    var_residual = anova_table{4,5};
    var_obj_ID = anova_table{3,5};
    obj_select_ind = var_obj_ID/(var_residual+var_obj_ID); % object selectivity index is between 0 and 1, the higher value, the more selective it is.
    Neuronlist(xxx).P_object_select_novel_end = anova_p(2);
    Neuronlist(xxx).obj_select_ind_novel_end = obj_select_ind;
    
    
    %% Recency - object selectivity interaction
    fractype4_5_ind = intersect(frac_sets.successfulfrac,find(ismember(Generaltask.Fractals(2,:), [4,5])));
    Generaltask.Fractals = IFI(Generaltask.Fractals, frac_sets); % calculate inter-fractal interval and save it back into the struct.
    recency_frac_sets = make_recency_frac_sets(Generaltask.Fractals(:,fractype4_5_ind));
    local_Factors = Generaltask.Fractals(:,fractype4_5_ind);
    local_FR = All_Fractal_FR(fractype4_5_ind);
    
    if numel(unique(local_Factors(1,recency_frac_sets.shortIFI_ind)))>=2 && numel(unique(local_Factors(1,recency_frac_sets.longIFI_ind)))>=2
        [~, anova_result]= regression_anova(local_FR(recency_frac_sets.shortIFI_ind,1), local_Factors(:,recency_frac_sets.shortIFI_ind), {'position','fractalID'});
        anova_p = anova_result.anova_p;
        anova_table = anova_result.anova_table;
        
        var_residual = anova_table{4,5};
        var_obj_ID = anova_table{3,5};
        obj_select_ind = var_obj_ID/(var_residual+var_obj_ID); % object selectivity index is between 0 and 1, the higher value, the more selective it is.
        Neuronlist(xxx).P_object_select_recent = anova_p(2);
        Neuronlist(xxx).obj_select_ind_recent = obj_select_ind;
        
        [~, anova_result]= regression_anova(local_FR(recency_frac_sets.longIFI_ind,1), local_Factors(:,recency_frac_sets.longIFI_ind), {'position','fractalID'});
        anova_p = anova_result.anova_p;
        anova_table = anova_result.anova_table;
        
        var_residual = anova_table{4,5};
        var_obj_ID = anova_table{3,5};
        obj_select_ind = var_obj_ID/(var_residual+var_obj_ID); % object selectivity index is between 0 and 1, the higher value, the more selective it is.
        Neuronlist(xxx).P_object_select_nonrecent = anova_p(2);
        Neuronlist(xxx).obj_select_ind_nonrecent = obj_select_ind;
        
    else
        Neuronlist(xxx).P_object_select_recent = nan;
        Neuronlist(xxx).obj_select_ind_recent = nan;
        Neuronlist(xxx).P_object_select_nonrecent = nan;
        Neuronlist(xxx).obj_select_ind_nonrecent = nan;
    end
    
    %% Sensory surprise - object selectivity interaction
    more_surprise_frac = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(2,:), [4,5]) & Generaltask.Fractals(3,:)==1));
    less_surprise_frac = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(2,:), [4,5]) & Generaltask.Fractals(3,:)==2));
    
    [~, anova_result]= regression_anova(All_Fractal_FR(more_surprise_frac,1), Generaltask.Fractals(:,more_surprise_frac), {'position','fractalID'});
    anova_p = anova_result.anova_p;
    anova_table = anova_result.anova_table;
    
    var_residual = anova_table{4,5};
    var_obj_ID = anova_table{3,5};
    obj_select_ind = var_obj_ID/(var_residual+var_obj_ID); % object selectivity index is between 0 and 1, the higher value, the more selective it is.
    Neuronlist(xxx).P_object_select_surprise = anova_p(2);
    Neuronlist(xxx).obj_select_ind_surprise = obj_select_ind;
    
    [~, anova_result]= regression_anova(All_Fractal_FR(less_surprise_frac,1), Generaltask.Fractals(:,less_surprise_frac), {'position','fractalID'});
    anova_p = anova_result.anova_p;
    anova_table = anova_result.anova_table;
    
    var_residual = anova_table{4,5};
    var_obj_ID = anova_table{3,5};
    obj_select_ind = var_obj_ID/(var_residual+var_obj_ID); % object selectivity index is between 0 and 1, the higher value, the more selective it is.
    Neuronlist(xxx).P_object_select_nonsurprise = anova_p(2);
    Neuronlist(xxx).obj_select_ind_nonsurprise = obj_select_ind;
    
end

end

function sq_err = Loss_exp_fun_cor(cor_matrix, exp_paras, binsize)
%%% calculate squared error given a set of parameters.
% cor_matrix here should be a lower trangle matrix, and each column should
% be the correlation coeffecient of a fixed time interval
exp_fun = @(k, A, B, t) A*(exp(-k*binsize/t)+B);
yvalue = exp_fun([1:size(cor_matrix,2)], exp_paras(1), exp_paras(2), exp_paras(3));

sqerr_matrix = (cor_matrix - yvalue).^2;
sq_err = nansum(sqerr_matrix(:));

end
