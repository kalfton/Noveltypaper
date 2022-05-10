%%%% find a good Novelty and recency encoding neuron and plot its SDF
example_neuron_outputfolder = fullfile(plotpath,'\sampleneurons');
raw_file_folder = '.\raw_data';
example_file_folder = '.\example_neuron_code\exampleneuron_raw_data';
Sampleneuron_set = find([Neuronlist_all.pred_nov_vs_fam]>0 & [Neuronlist_all.P_pred_nov_vs_fam]<StatisticalThreshold &...
    [Neuronlist_all.pred_vs_unpred_fam]>0 & [Neuronlist_all.P_pred_vs_unpred_fam_perm]<StatisticalThreshold &...
    [Neuronlist_all.recency_ind_match_pos]>0 & [Neuronlist_all.P_recency_ind_match_pos]<StatisticalThreshold); % violation_ind_4roc>0 & Pviolation_ind_4roc<StatisticalThreshold

yaxislim = [-1,3];
exampleneurons =[7,18];

plot_nrow = 5;
plot_ncol= 5;
plotplacesetx = {1,2,3,4,5};
plotplacesety = {1,1,1,1,1};

title_for_plot = {'Novelty','Sensory surprise', 'Recency'};
SDF_for_plot = { 'pred_nov_vs_fam_sdfs', 'pred_vs_unpred_fam_sdfs', 'recency_ind_match_pos_sdfs'};
indices_in_struct = {'pred_nov_vs_fam', 'pred_vs_unpred_fam', 'recency_ind_match_pos'};
Pvalues_in_struct = {'P_pred_nov_vs_fam', 'P_pred_vs_unpred_fam_perm', 'P_recency_ind_match_pos'};

for i = 1:numel(exampleneurons)
    %% calculate the novelty index, recency index and surprise index for one neuron:
    id = exampleneurons(i);
    sampleneuron = Neuronlist_all(Sampleneuron_set(id));
    
    clear(sampleneuron.name, 'Generaltask');
    %load the raw neuron file
    try
        if strcmpi(sampleneuron.monkeyName, 'S')
            load(fullfile(raw_file_folder, 'Monkey_S_raw', sampleneuron.filename), sampleneuron.name, 'Generaltask');
        elseif strcmpi(sampleneuron.monkeyName, 'L')
            load(fullfile(raw_file_folder, 'Monkey_L_raw', sampleneuron.filename), sampleneuron.name, 'Generaltask');
        end
    catch % use the local example file, if the full raw file is not available
        load(fullfile(example_file_folder, sampleneuron.filename), sampleneuron.name, 'Generaltask');
    end
    
    baselinesubtraction = 'None';
    clear currentstruct;
    
    channelname = sampleneuron.name;
    currentstruct.name = sampleneuron.name;
    currentstruct.filename = sampleneuron.filename;
    %%
    rng(0);
    %shuffling_num = 10000;
    Smoothing=3; %1 - box car; 2- gaus; 3 - epsp (causal); specs of kernel are defined in Smooth_Histogram function
    gauswindow_ms=50;
    Fractalsdfwindow = 200:800;
    
    trialsplitlogical = true(size(Generaltask.successtrial));
    
    fracsplitlogical = kron(trialsplitlogical',true(1,3));
    
    successtrial_split = find(Generaltask.successtrial == 1 & trialsplitlogical);
    
    %% generate fractal sets, put it outside the for loop.
    frac_sets = make_frac_sets(Generaltask, fracsplitlogical);
    
    zeropoint_f = 300;
    fractaldur = 1:500;
    Hist_raster_f = [];
    
    %% get rid of repeated 7999 fractal
    for i = 1:size(Generaltask.Fractals,2)
        if(Generaltask.Fractals(1,i)==7999)
            switch Generaltask.Fractals(3,i)
                case 2
                    if(ismember(Generaltask.Fractals(1,i-1),[-1,7999]))
                        Generaltask.Fractals(1,i)=-1;
                    end
                case 3
                    if(ismember(Generaltask.Fractals(1,i-1),[-1,7999]) || ismember(Generaltask.Fractals(1,i-2),[-1,7999]))
                        Generaltask.Fractals(1,i)=-1;
                    end
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eval(['fracsptimes = ' channelname '.fracsptimes;'])
    for x=1:size(Generaltask.Fractals,2)
        Hist_template_f = -zeropoint_f:(500+100);%surpass the span a little bit
        Hist_raw_f = histcounts((fracsptimes{1,x}(1:end)-Generaltask.Fractals(4,x))*1000, Hist_template_f);
        Hist_raster_f = [Hist_raster_f; Hist_raw_f];
    end
    
    Hist_sdf_f=plot_mean_psth({Hist_raster_f},gauswindow_ms,1,size(Hist_raster_f,2),1);
    Rasters_f = Hist_raster_f;
    Rasterscs_f=Hist_sdf_f;
    Hist_template_f = Hist_template_f(1:end-1);
    
    ObjectSpikeFR = nansum(  Rasters_f( :, zeropoint_f+fractaldur ),2)'*1000/length(fractaldur) ;
    PreObject=nansum(Rasters_f( :, 220:320 ),2)'*1000/101;
    %Fix_aveFR = Fixation_trtypemean(Generaltask.Fractals(2,:));
    if strcmpi(baselinesubtraction,'None')
        All_Fractal_FR=ObjectSpikeFR; All_Fractal_FR=All_Fractal_FR';
    else
        warning('wrong variable of baselinesubtraction variable, set it to default ''None'' ');
        All_Fractal_FR=ObjectSpikeFR; All_Fractal_FR=All_Fractal_FR';
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in type 4 and 5 recency effect
    
    Generaltask.Fractals = IFI(Generaltask.Fractals, frac_sets); % calculate inter-fractal interval and save it back into the struct.
    
    fractype4_5_ind = intersect(frac_sets.successfulfrac,find(ismember(Generaltask.Fractals(2,:), [4,5])));
    Alltimedifference = Generaltask.Fractals(6,fractype4_5_ind);
    [Alltimedifference, I] = sort(Alltimedifference);
    fractype4_5_ind = fractype4_5_ind(I);
    
    % Get rid of the first time appearance
    first_appearance4_5_ind = fractype4_5_ind(Alltimedifference==inf);
    fractype4_5_ind = fractype4_5_ind(Alltimedifference~=inf);
    Alltimedifference = Alltimedifference(Alltimedifference~=inf);
    
    % The trials are splited for cross validation
    [fractype4_5_ind, I] = intersect(fractype4_5_ind, frac_sets.successfulfrac_split);
    Alltimedifference = Alltimedifference(I);
    [Alltimedifference, I] = sort(Alltimedifference);
    fractype4_5_ind = fractype4_5_ind(I);
    
    
    %Get corresponding firing rate
    FiringR_45 = All_Fractal_FR(fractype4_5_ind,1);
    
    
    FiringR_45_subtracted = regression_anova(FiringR_45, Generaltask.Fractals(:,fractype4_5_ind), {'position', 'fractalID'});
    
    
    %% recency index
    recency_frac_sets = make_recency_frac_sets(Generaltask.Fractals(:,fractype4_5_ind));
    Rasterscs_f_=Rasterscs_f(fractype4_5_ind,:);
    
    if ~isempty(recency_frac_sets.shortIFI_ind)
        currentstruct.recency_ind_match_pos = rocarea3(FiringR_45_subtracted(recency_frac_sets.shortIFI_ind),FiringR_45_subtracted(recency_frac_sets.longIFI_ind));
        currentstruct.P_recency_ind_match_pos = ranksum(FiringR_45_subtracted(recency_frac_sets.shortIFI_ind),FiringR_45_subtracted(recency_frac_sets.longIFI_ind));
    else
        currentstruct.recency_ind_match_pos = nan;
        currentstruct.P_recency_ind_match_pos = nan;
    end
    currentstruct.recency_ind_match_pos_sdfs =[nanmean(Rasterscs_f_(recency_frac_sets.shortIFI_ind,Fractalsdfwindow));
        nanmean(Rasterscs_f_(recency_frac_sets.longIFI_ind,Fractalsdfwindow))];
    
    
    %bootstrapping to get error bars
    if ~isempty(recency_frac_sets.shortIFI_ind)
        FiringR_shortIFI = FiringR_45_subtracted(recency_frac_sets.shortIFI_ind);
        FiringR_longIFI = FiringR_45_subtracted(recency_frac_sets.longIFI_ind);
        recency_ind_match_pos_shuffled = zeros(shuffling_num,1);
        for ii = 1:shuffling_num
            recency_ind_match_pos_shuffled(ii) = rocarea3(datasample(FiringR_shortIFI,numel(FiringR_shortIFI)), datasample(FiringR_longIFI,numel(FiringR_longIFI)));
        end
        recency_ind_match_pos_shuffled = sort(recency_ind_match_pos_shuffled);
        currentstruct.recency_ind_match_pos_upperbound = recency_ind_match_pos_shuffled(ceil(shuffling_num - shuffling_num*0.005));
        currentstruct.recency_ind_match_pos_lowerbound = recency_ind_match_pos_shuffled(floor(shuffling_num*0.005));
        currentstruct.recency_ind_match_pos_std = std(recency_ind_match_pos_shuffled);
    else
        currentstruct.recency_ind_match_pos = nan;
        currentstruct.P_recency_ind_match_pos = nan;
        currentstruct.recency_ind_match_pos_std = nan;
    end
    
    
    
    % predicted novel vs predicted familiar
    fractrials = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6401,6404,7999])));
    fractrials = intersect(fractrials, find(ismember(Generaltask.Fractals(2,:), [2,3,6])));
    fractrials = intersect(fractrials, find(ismember(Generaltask.Fractals(3,:), 2)));
    FiringR = All_Fractal_FR(fractrials,1);
    
    FR_1 = FiringR(Generaltask.Fractals(2,fractrials)==6);
    FR_2 = FiringR(ismember(Generaltask.Fractals(2,fractrials),[2,3]));
    
    currentstruct.pred_nov_vs_fam = rocarea3(FR_2,FR_1);
    currentstruct.P_pred_nov_vs_fam = ranksum(FiringR(ismember(Generaltask.Fractals(2,fractrials),[2,3])),FiringR(Generaltask.Fractals(2,fractrials)==6));
    Rasterscs_f_=Rasterscs_f(fractrials,:);
    currentstruct.pred_nov_vs_fam_sdfs =[nanmean(Rasterscs_f_(ismember(Generaltask.Fractals(2,fractrials),[2,3]),Fractalsdfwindow));
        nanmean(Rasterscs_f_(Generaltask.Fractals(2,fractrials)==6,Fractalsdfwindow))];
    
    %bootstrapping to get error bars
    FiringR_PF = FR_2;
    FiringR_PN = FR_1;
    pred_nov_vs_fam_shuffled = zeros(shuffling_num,1);
    for ii = 1:shuffling_num
        pred_nov_vs_fam_shuffled(ii) = rocarea3(datasample(FiringR_PF,numel(FiringR_PF)), datasample(FiringR_PN,numel(FiringR_PN)));
    end
    pred_nov_vs_fam_shuffled = sort(pred_nov_vs_fam_shuffled);
    currentstruct.pred_nov_vs_fam_upperbound = pred_nov_vs_fam_shuffled(ceil(shuffling_num - shuffling_num*0.005));
    currentstruct.pred_nov_vs_fam_lowerbound = pred_nov_vs_fam_shuffled(floor(shuffling_num*0.005));
    currentstruct.pred_nov_vs_fam_std = std(pred_nov_vs_fam_shuffled);
    
    
    %% calculating surprise index
    % prepare the firing rate to be recency adapted.
    familiar45type = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6500:6510])));
    FR_45 = All_Fractal_FR(familiar45type,1);
    [FR45_subtracted, anova_result]= regression_anova(FR_45, Generaltask.Fractals(:,familiar45type), {'recency'});
    Al_FR_subtracted = All_Fractal_FR;
    Al_FR_subtracted(familiar45type,1) = FR45_subtracted;
    
    fractrials = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6400:6405, 6500:6505])));
    fractrials = intersect(fractrials, find(ismember(Generaltask.Fractals(2,:), [2,3,4,5]))); % trial type
    %FiringR = All_Fractal_FR(fractrials,1);
    FiringR = Al_FR_subtracted(fractrials,1);
    Factor_x = zeros(length(fractrials), 2);
    Factor_x(:,1) = Generaltask.Fractals(3,fractrials)'; % order position
    Factor_x(:,2) = Generaltask.Fractals(1,fractrials)'; % fractal ID
    
    % Average ROC's
    FiringR_3rd_23 = FiringR(Factor_x(:,1) == 3 & ismember(Factor_x(:,2),6400:6405));
    FiringR_3rd_45 = FiringR(Factor_x(:,1) == 3 & ismember(Factor_x(:,2),6500:6505));
    currentstruct.pred_vs_unpred_fam = rocarea3(FiringR_3rd_23, FiringR_3rd_45);
    currentstruct.P_pred_vs_unpred_fam_perm = ranksum(FiringR_3rd_23, FiringR_3rd_45);
    %
    Rasterscs_f_=Rasterscs_f(fractrials,:);
    currentstruct.pred_vs_unpred_fam_sdfs =[ nanmean(Rasterscs_f_(Factor_x(:,1) == 3 & ismember(Factor_x(:,2),6400:6405),Fractalsdfwindow))
        nanmean(Rasterscs_f_(Factor_x(:,1) == 3 & ismember(Factor_x(:,2),6500:6505),Fractalsdfwindow))    ];
    
    %bootstrapping to get error bars
    pred_vs_unpred_fam_shuffled = zeros(shuffling_num,1);
    for ii = 1:shuffling_num
        pred_vs_unpred_fam_shuffled(ii) = rocarea3(datasample(FiringR_3rd_23,numel(FiringR_3rd_23)), datasample(FiringR_3rd_45,numel(FiringR_3rd_45)));
    end
    pred_vs_unpred_fam_shuffled = sort(pred_vs_unpred_fam_shuffled);
    currentstruct.pred_vs_unpred_fam_upperbound = pred_vs_unpred_fam_shuffled(ceil(shuffling_num - shuffling_num*0.005));
    currentstruct.pred_vs_unpred_fam_lowerbound = pred_vs_unpred_fam_shuffled(floor(shuffling_num*0.005));
    currentstruct.pred_vs_unpred_fam_std = std(pred_vs_unpred_fam_shuffled);
    
    %normalize the indices
    currentstruct.pred_nov_vs_fam = 2*(currentstruct.pred_nov_vs_fam-0.5);
    currentstruct.pred_vs_unpred_fam = 2*(currentstruct.pred_vs_unpred_fam-0.5);
    currentstruct.recency_ind_match_pos = 2*(currentstruct.recency_ind_match_pos-0.5);
    
    
    % end of calculating the indices
    %% Plot the example neuron SDF
    figure;
    for xy = 1:length(title_for_plot)
        
        nsubplot(plot_nrow, plot_ncol, plotplacesety{xy}, plotplacesetx{xy});
        
        title([title_for_plot{xy}]);
        plot(currentstruct.(SDF_for_plot{xy})(1,:),'-b', 'LineWidth', 1);
        plot(currentstruct.(SDF_for_plot{xy})(2,:),'-r', 'LineWidth', 1);
        if xy==1
            localyaxislim = ylim;
        end
        
        text(400,0, ['p=' mat2str(currentstruct.(Pvalues_in_struct{xy}),3)]);
        text(400,20, ['ind=' mat2str(currentstruct.(indices_in_struct{xy}),3)]);
        
        xlim([0,600]);
        ylim(localyaxislim);
        y_lim = get(gca,'ylim');
        plot([100,100],y_lim,'k');
    end
    % electrode ID and filename
    xy = 5;
    nsubplot(plot_nrow, plot_ncol, plotplacesety{xy}, plotplacesetx{xy});
    text(0,1,currentstruct.filename, 'interpreter', 'none');
    text(0,2,currentstruct.name);
    ylim([0,4]);
    axis off;
    
    
    %% plot the bar plot
    xy=4;
    nsubplot(plot_nrow, plot_ncol, plotplacesety{xy}, plotplacesetx{xy});
    
    bar([1,2,3], [currentstruct.pred_nov_vs_fam, currentstruct.pred_vs_unpred_fam, currentstruct.recency_ind_match_pos], 'r');
    
    errorbar([1,2,3], [currentstruct.pred_nov_vs_fam, currentstruct.pred_vs_unpred_fam, currentstruct.recency_ind_match_pos], ...
        [currentstruct.pred_nov_vs_fam_std, currentstruct.pred_vs_unpred_fam_std, currentstruct.recency_ind_match_pos_std],...
        '.');
    
    set(gca, 'xtick', [1,2,3], 'xticklabel', {'Novelty', 'Sensory surprise', 'Recency'});
    xlim([0.5,3.5]);
    ylabel('Index value');
    ylim([0,1]);
    
    set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
    print(gcf,'-dpdf', '-painters',fullfile(example_neuron_outputfolder, ['sample neuoron' sampleneuron.name, sampleneuron.filename '.pdf']));
    
end
