%%%% find a good Novelty and recency encoding neuron and plot its SDF
sample_neuron_folder = '.\plots\sampleneurons';
Sampleneuron_set = find([Neuronlist_good.pred_nov_vs_fam]>0 & [Neuronlist_good.P_pred_nov_vs_fam]<StatisticalThreshold &...
     [Neuronlist_good.pred_vs_unpred_fam]>0 & [Neuronlist_good.P_pred_vs_unpred_fam_perm]<StatisticalThreshold &...
     [Neuronlist_good.recency_ind_match_pos]>0 & [Neuronlist_good.P_recency_ind_match_pos]<StatisticalThreshold); % violation_ind_4roc>0 & Pviolation_ind_4roc<StatisticalThreshold

yaxislim = [-1,3];
id =7;
%sampleneuron = Neuronlist_good(Sampleneuron_set(id));
%Neurondata = vertcat(Neuronlist_good(signfid).pred_nov_vs_fam_sdfs);
% figure; 
% 
% subplot(2,2,1)
% hold on;
% plot(sampleneuron.pred_nov_vs_fam_sdfs(1,:),'-b', 'LineWidth', 1);
% plot(sampleneuron.pred_nov_vs_fam_sdfs(2,:),'-r', 'LineWidth', 1);
% xlim([0,600]);
% ylim(yaxislim);
% y_lim = get(gca,'ylim');
% plot([100,100],y_lim,'k');
% title('Novelty');
% 
% subplot(2,2,2)
% hold on;
% plot(sampleneuron.recency_ind_sdfs(1,:),'-b', 'LineWidth', 1);
% plot(sampleneuron.recency_ind_sdfs(2,:),'-r', 'LineWidth', 1);
% title('Recency');


%% automatically switch lines.

title_for_plot = {'Novelty', 'Recency','Sensory surprise', 'Violation'};
SDF_for_plot = { 'pred_nov_vs_fam_sdfs', 'recency_ind__match_pos_sdfs', 'pred_vs_unpred_fam_sdfs', 'violation_ind_comp_sdf'};
indices_in_struct = {'pred_nov_vs_fam', 'recency_ind_match_pos', 'pred_vs_unpred_fam', 'violation_ind_4roc'};
Pvalues_in_struct = {'P_pred_nov_vs_fam', 'P_recency_ind_match_pos', 'P_pred_vs_unpred_fam_perm', 'P_violation_ind_4roc'};
plotplacesetx = {10:30, 40:60,70:90,100:120,130:150};
plotplacesety = {5:25, 5:25, 5:25, 5:25, 5:25};
ygap = 30

figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
figurenum  = 1;
rownum = 0;
for xx = 1:length(Sampleneuron_set)
    if (rownum+1)*ygap>169
        % save current figure
        print(gcf,'-dpdf', '-painters',fullfile(sample_neuron_folder, ['sampleneurons' mat2str(figurenum) '.pdf']));
        % new figure
        figurenum = figurenum+1;
        figure;
        set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
        rownum = 0;
    end
    for xy = 1:length(title_for_plot)
        
        sampleneuron = Neuronlist_good(Sampleneuron_set(xx));
        
        nsubplot(169,169, rownum*ygap+plotplacesety{xy}, plotplacesetx{xy});
        
        title([title_for_plot{xy}]);
        datastd = sampleneuron.zscorestd;
        datamean = sampleneuron.zscoremean;
        plot(datamean+datastd*sampleneuron.(SDF_for_plot{xy})(1,:),'-b', 'LineWidth', 1);
        plot(datamean+datastd*sampleneuron.(SDF_for_plot{xy})(2,:),'-r', 'LineWidth', 1);
        if xy==1
            localyaxislim = yaxislim*datastd+datamean;
        end
        
        text(400,2, ['p=' mat2str(sampleneuron.(Pvalues_in_struct{xy}),3)]);
        text(400,2.5, ['ind=' mat2str(sampleneuron.(indices_in_struct{xy}),3)]);
        
        xlim([0,600]);
        ylim(localyaxislim);
        y_lim = get(gca,'ylim');
        plot([100,100],y_lim,'k');
    end
    % electrode ID and filename
    xy = 5;
    nsubplot(169,169, rownum*ygap+plotplacesety{xy}, plotplacesetx{xy});
    text(0,1,sampleneuron.filename);
    text(0,2,['Depth ' mat2str(sampleneuron.depth,3)]);
    text(0,3,sampleneuron.name);
    ylim([0,4]);
    axis off;
    
    rownum = rownum+1;
end
print(gcf,'-dpdf', '-painters',fullfile([sample_neuron_folder, 'sampleneurons' mat2str(figurenum) '.pdf']));

%%
% plot the bar plot of recency index and surprise index for one neuron
sampleneuron = Neuronlist_good(Sampleneuron_set(id));

% error bars should be caluculated by bootstrapping.
% load the raw neuron file
load(fullfile('Y:\PLEXON_GRAYARRAY_Slayer\NovelFractalLearning_Slayer_combinedsorting_temp', sampleneuron.filename), sampleneuron.name, 'Generaltask');
zscorestd = sampleneuron.zscorestd;
zscoremean = sampleneuron.zscoremean;
channelname = sampleneuron.name;

baselinesubtraction = 'None';
clear currentstruct;


% copy the alg. from CreateIndices_V3
%%
rng(0);
%shuffling_num = 10000;
Smoothing=3; %1 - box car; 2- gaus; 3 - epsp (causal); specs of kernel are defined in Smooth_Histogram function
gauswindow_ms=50;

trialsplitlogical = true(size(Generaltask.successtrial));

fracsplitlogical = kron(trialsplitlogical',true(1,3));

successtrial_split = find(Generaltask.successtrial == 1 & trialsplitlogical);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% successful fractal & violation fractal in type 2&3
% in trial type 2&3, pick the same number of normal fractal as
% comparison
successfulfrac = find(~isnan(Generaltask.Fractals(5,:)));
successfulfrac_split = find(~isnan(Generaltask.Fractals(5,:))& fracsplitlogical);
numberb=2;

temp_ind = intersect(successfulfrac_split,find(ismember(Generaltask.Fractals(1,:),[6401])));
violation_frac_2nd_6401 = temp_ind(Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1);
temp_ind = setdiff(temp_ind,violation_frac_2nd_6401);
normal_frac_2nd_contrast_6401 = temp_ind(randperm(numel(temp_ind)));
normal_frac_2nd_contrast_6401 = sort(normal_frac_2nd_contrast_6401(1:min(numberb*numel(violation_frac_2nd_6401),numel(normal_frac_2nd_contrast_6401))));

temp_ind = intersect(successfulfrac_split,find(ismember(Generaltask.Fractals(1,:),[6402])));
violation_frac_3rd_6402 = temp_ind((Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1) & (Generaltask.Fractals(1,temp_ind-2) ~= Generaltask.Fractals(1,temp_ind)-2)); % 3nd violation fractal
violation_frac_3rd_2_6402 = temp_ind(Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1 & Generaltask.Fractals(1,temp_ind-2) == Generaltask.Fractals(1,temp_ind)-2); % 3rd fractal following the 2nd violation fractal
temp_ind = setdiff(temp_ind,violation_frac_3rd_6402);
normal_frac_3rd_contrast_6402 = temp_ind(randperm(numel(temp_ind)));
try
    normal_frac_3rd_contrast_6402 = sort(normal_frac_3rd_contrast_6402(1:numberb*numel(violation_frac_3rd_6402)));
catch %some session happened to have a lot of violated fractals and the normal fractal is not enough to match them, here is the solution
    normal_frac_3rd_contrast_6402 = sort(normal_frac_3rd_contrast_6402(1:floor(length(normal_frac_3rd_contrast_6402))));
end

%%%%%%%%%%%%%%%%%%%%

temp_ind = intersect(successfulfrac_split,find(ismember(Generaltask.Fractals(1,:),[6404])));
violation_frac_2nd_6404 = temp_ind(Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1);
temp_ind = setdiff(temp_ind,violation_frac_2nd_6404);
normal_frac_2nd_contrast_6404 = temp_ind(randperm(numel(temp_ind)));
normal_frac_2nd_contrast_6404 = sort(normal_frac_2nd_contrast_6404(1:numberb*numel(violation_frac_2nd_6404)));

temp_ind = intersect(successfulfrac_split,find(ismember(Generaltask.Fractals(1,:),[6405])));
violation_frac_3rd_6405 = temp_ind((Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1) & (Generaltask.Fractals(1,temp_ind-2) ~= Generaltask.Fractals(1,temp_ind)-2)); % 3nd violation fractal
violation_frac_3rd_2_6405 = temp_ind(Generaltask.Fractals(1,temp_ind-1) ~= Generaltask.Fractals(1,temp_ind)-1 & Generaltask.Fractals(1,temp_ind-2) == Generaltask.Fractals(1,temp_ind)-2); % 3rd fractal following the 2nd violation fractal
temp_ind = setdiff(temp_ind,violation_frac_3rd_6405);
normal_frac_3rd_contrast_6405 = temp_ind(randperm(numel(temp_ind)));
normal_frac_3rd_contrast_6405 = sort(normal_frac_3rd_contrast_6405(1:numberb*numel(violation_frac_3rd_6405)));

violation_frac_2nd = [violation_frac_2nd_6401,violation_frac_2nd_6404];
violation_frac_3rd = [violation_frac_3rd_6402,violation_frac_3rd_6405];
violation_frac_3rd_2 = [violation_frac_3rd_2_6402,violation_frac_3rd_2_6405];
normal_frac_2nd_contrast = [normal_frac_2nd_contrast_6401,normal_frac_2nd_contrast_6404];
normal_frac_3rd_contrast = [normal_frac_3rd_contrast_6402,normal_frac_3rd_contrast_6405];

%%%%%%%%%%%%%%%%%%%%

%violation_frac = intersect([violation_frac_2nd, violation_frac_3rd],successfulfrac);
successfulfrac_noviol = setdiff(successfulfrac_split,[violation_frac_2nd, violation_frac_3rd, violation_frac_3rd_2,normal_frac_3rd_contrast,normal_frac_2nd_contrast],'sorted');

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
Hist_sdf_f = (Hist_sdf_f-zscoremean)./zscorestd;
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in type 4 and 5 recency effect
Fractalsdfwindow = 200:800;

% Time intervals between fractals may also affect neuron's firing rate.
fractype4_5 = cell(6,1);
timedifference = cell(6,1);
for i=1:6
    fractype4_5{i} = sort(intersect(successfulfrac,find(Generaltask.Fractals(1,:)==(6500+i-1)))); % from 6500 to 6505
    fracnum = fractype4_5{i};
    timedifference{i} = Generaltask.Fractals(4,fracnum(2:end))-Generaltask.Fractals(4,fracnum(1:end-1)); % caculate time gap between the appearance of fractal
    timedifference{i} = [inf, timedifference{i}]; % the first time the fractal appears, set the time gap to be inf
end

[Alltimedifference, I] = sort(cell2mat(timedifference'));
fractype4_5_ind = cell2mat(fractype4_5');
fractype4_5_ind = fractype4_5_ind(I);

% Get rid of the first time appearance
first_appearance4_5_ind = fractype4_5_ind(Alltimedifference==inf);
fractype4_5_ind = fractype4_5_ind(Alltimedifference~=inf);
Alltimedifference = Alltimedifference(Alltimedifference~=inf);

% The trials are splited for cross validation
[fractype4_5_ind, I] = intersect(fractype4_5_ind, successfulfrac_split);
Alltimedifference = Alltimedifference(I);
[Alltimedifference, I] = sort(Alltimedifference);
fractype4_5_ind = fractype4_5_ind(I);


%Get corresponding firing rate
FiringR_45 = All_Fractal_FR(fractype4_5_ind,1);

% separate the time difference into n>=3 bins:
% Fractal_timedifference store the time difference, the column
% dimension equals all fractals,  the first column stores time
% difference and the second stores category: short(1st bin) medium(2nd
% to n-1 bins) long(last bin)
binnum = 5; % you can change the bin number here
timedifference_category = zeros(size(Alltimedifference));
timedifference_category(1:floor(length(timedifference_category)/binnum)) = 1;
timedifference_category((floor(length(timedifference_category)/binnum)+1):(floor(length(timedifference_category)*(binnum-1)/binnum))) = 2;
timedifference_category((floor(length(timedifference_category)*(binnum-1)/binnum)+1):end) = 3;

time_crit_1 = Alltimedifference(floor(length(timedifference_category)/binnum));
time_crit_2 = Alltimedifference(floor(length(timedifference_category)*(binnum-1)/binnum));
time_crit_3 = Alltimedifference(end);

Factor_x_45 = zeros(length(fractype4_5_ind), 3);
Factor_x_45(:,1) = Generaltask.Fractals(3,fractype4_5_ind)'; % order factor
Factor_x_45(:,2) = Generaltask.Fractals(1,fractype4_5_ind)'; % fractal ID factor
Factor_x_45(:,3) = timedifference_category';


[anova_p,anova_table,anova_stats] = anovan(FiringR_45, Factor_x_45(:,1:2),'display','off','model','linear','varnames',{'position','fractalID'});
positioneff = anova_stats.coeffs(contains(anova_stats.coeffnames, {'position'}));% 3 position
fractalIDeff = anova_stats.coeffs(contains(anova_stats.coeffnames, {'fractalID'}));% 6 different fractals
%recencyeff = anova_stats.coeffs(contains(anova_stats.coeffnames, {'recency'}));

%fractalIDeff should have 6 items. if it is less than we should remap
%it
if length(fractalIDeff)<6
    uniquefractalID = sort(unique(Factor_x_45(~isnan(FiringR_45),2)));
    oldfractalIDeff = fractalIDeff;
    fractalIDeff=zeros(6,1);
    fractalIDeff(uniquefractalID-6500+1)=oldfractalIDeff;
end
% recency index & p value
FiringR_45_subtracted = FiringR_45-positioneff(Generaltask.Fractals(3,fractype4_5_ind))-fractalIDeff(Generaltask.Fractals(1,fractype4_5_ind)-6500+1); % Firing rate with position effect and fractal ID effect subtracted

currentstruct.recency_ind = rocarea3(FiringR_45_subtracted(Factor_x_45(:,3)==1),FiringR_45_subtracted(Factor_x_45(:,3)==3));
currentstruct.P_recency_ind = ranksum(FiringR_45_subtracted(Factor_x_45(:,3)==1),FiringR_45_subtracted(Factor_x_45(:,3)==3));

%     rocarea3(FiringR_45(Factor_x_45(:,3)==1),FiringR_45(Factor_x_45(:,3)==3))
%     ranksum(FiringR_45(Factor_x_45(:,3)==1),FiringR_45(Factor_x_45(:,3)==3))


Rasterscs_f_=Rasterscs_f(fractype4_5_ind,:);
currentstruct.recency_ind_sdfs =[nanmean(Rasterscs_f_(Factor_x_45(:,3)==1,Fractalsdfwindow));
    nanmean(Rasterscs_f_(Factor_x_45(:,3)==3,Fractalsdfwindow))];
currentstruct.P_recency_ind_shuffle = currentstruct.P_recency_ind; % shuffling of AUC is the same as ranksum test
%% second way to calculate recency index
% Short IFI fractal: fractals repeating inside trial
% Long IFI fractal: fractals repeating outside trial
% (matching position)
shortIFI_ind_2nd = find(Alltimedifference<1.5 & Factor_x_45(:,1)'==2);
shortIFI_ind_3rd = find(Alltimedifference<1.5 & Factor_x_45(:,1)'==3);
num_short_2nd_pos = numel(shortIFI_ind_2nd);
num_short_3rd_pos = numel(shortIFI_ind_3rd);

longIFI_ind_2nd = find(Alltimedifference>1.5 & Factor_x_45(:,1)'==2);
longIFI_ind_3rd = find(Alltimedifference>1.5 & Factor_x_45(:,1)'==3);
num_long_2nd_pos = numel(longIFI_ind_2nd);
num_long_3rd_pos = numel(longIFI_ind_3rd);

num_2nd_pos = min([num_short_2nd_pos, num_long_2nd_pos, num_short_3rd_pos, num_long_3rd_pos]);
num_3rd_pos = num_2nd_pos;

shortIFI_ind_2nd = shortIFI_ind_2nd(1:num_2nd_pos);
shortIFI_ind_3rd = shortIFI_ind_3rd(1:num_3rd_pos);
longIFI_ind_2nd = longIFI_ind_2nd(end-num_2nd_pos+1:end);
longIFI_ind_3rd = longIFI_ind_3rd(end-num_3rd_pos+1:end);

shortIFI_ind = [shortIFI_ind_2nd, shortIFI_ind_3rd];
longIFI_ind = [longIFI_ind_2nd, longIFI_ind_3rd];

if ~isempty(shortIFI_ind)
    currentstruct.recency_ind_match_pos = rocarea3(FiringR_45_subtracted(shortIFI_ind),FiringR_45_subtracted(longIFI_ind));
    currentstruct.P_recency_ind_match_pos = ranksum(FiringR_45_subtracted(shortIFI_ind),FiringR_45_subtracted(longIFI_ind));
    
    %bootstrapping to get error bars
    FiringR_shortIFI = FiringR_45_subtracted(shortIFI_ind);
    FiringR_longIFI = FiringR_45_subtracted(longIFI_ind);
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

currentstruct.recency_ind__match_pos_sdfs =[nanmean(Rasterscs_f_(shortIFI_ind,Fractalsdfwindow));
    nanmean(Rasterscs_f_(longIFI_ind,Fractalsdfwindow))];



% predicted novel vs predicted familiar
fractrials = intersect(successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6401,6404,7999])));
fractrials = intersect(fractrials, find(ismember(Generaltask.Fractals(2,:), [2,3,6])));
fractrials = intersect(fractrials, find(ismember(Generaltask.Fractals(3,:), 2)));
FiringR = All_Fractal_FR(fractrials,1);
FR_1 = FiringR(Generaltask.Fractals(2,fractrials)==6);
FR_2 = FiringR(ismember(Generaltask.Fractals(2,fractrials),[2,3]));

try
    currentstruct.pred_nov_vs_fam = rocarea3(FiringR(ismember(Generaltask.Fractals(2,fractrials),[2,3])),FiringR(Generaltask.Fractals(2,fractrials)==6));
    currentstruct.P_pred_nov_vs_fam = ranksum(FiringR(ismember(Generaltask.Fractals(2,fractrials),[2,3])),FiringR(Generaltask.Fractals(2,fractrials)==6));
    
    Rasterscs_f_=Rasterscs_f(fractrials,:);
    currentstruct.pred_nov_vs_fam_sdfs =[nanmean(Rasterscs_f_(ismember(Generaltask.Fractals(2,fractrials),[2,3]),Fractalsdfwindow));
        nanmean(Rasterscs_f_(Generaltask.Fractals(2,fractrials)==6,Fractalsdfwindow))];
catch
    currentstruct.pred_nov_vs_fam = nan;
    currentstruct.P_pred_nov_vs_fam = nan;
    currentstruct.pred_nov_vs_fam_sdfs = nan;
end

%bootstrapping to get error bars
FiringR_PF = FiringR(ismember(Generaltask.Fractals(2,fractrials),[2,3]));
FiringR_PN = FiringR(Generaltask.Fractals(2,fractrials)==6);
pred_nov_vs_fam_shuffled = zeros(shuffling_num,1);
for ii = 1:shuffling_num
    pred_nov_vs_fam_shuffled(ii) = rocarea3(datasample(FiringR_PF,numel(FiringR_PF)), datasample(FiringR_PN,numel(FiringR_PN)));
end
pred_nov_vs_fam_shuffled = sort(pred_nov_vs_fam_shuffled);
currentstruct.pred_nov_vs_fam_upperbound = pred_nov_vs_fam_shuffled(ceil(shuffling_num - shuffling_num*0.005));
currentstruct.pred_nov_vs_fam_lowerbound = pred_nov_vs_fam_shuffled(floor(shuffling_num*0.005));
currentstruct.pred_nov_vs_fam_std = std(pred_nov_vs_fam_shuffled);


% predicted familar vs unpredicted familiar
fractrials = intersect(successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6400:6405, 6500:6505])));
fractrials = intersect(fractrials, find(ismember(Generaltask.Fractals(2,:), [2,3,4,5]))); % trial type
FiringR = All_Fractal_FR(fractrials,1);
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

%% new way of calculating surprise index
%% regress out recency
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



%% end of copying the code

figure;
xy=1;
nsubplot(169,169, ygap+plotplacesety{xy}, plotplacesetx{xy});

bar([1,2,3], 2*[currentstruct.pred_nov_vs_fam-0.5, currentstruct.pred_vs_unpred_fam-0.5, currentstruct.recency_ind_match_pos-0.5], 'r');

errorbar([1,2,3], 2*[currentstruct.pred_nov_vs_fam-0.5, currentstruct.pred_vs_unpred_fam-0.5, currentstruct.recency_ind_match_pos-0.5], ...
    [currentstruct.pred_nov_vs_fam_std, currentstruct.pred_vs_unpred_fam_std, currentstruct.recency_ind_match_pos_std],...
     '.');

set(gca, 'xtick', [1,2,3], 'xticklabel', {'Novelty', 'Sensory surprise', 'Recency'});
xlim([0.5,3.5]);
ylabel('Index value');
ylim([0,0.8]);

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(sample_neuron_folder, ['sample neuoron bar plot.pdf']));


%% population averaged sdfs
% Sampleneuron_set = find(pred_nov_vs_fam>0 & Ppred_nov_vs_fam'<StatisticalThreshold &...
%     pred_vs_unpred_fam>0 & Ppred_vs_unpred_fam'<StatisticalThreshold &...
%     recency_ind>0 & Precency_ind'<StatisticalThreshold);

Sampleneuron_set = find([Neuronlist_good.pred_nov_vs_fam]>0 & [Neuronlist_good.P_pred_nov_vs_fam]<StatisticalThreshold);
yaxislim = [-0.5,1.0];

figure;
for xy = 1:length(title_for_plot)
    
    %sampleneuron = Neuronlist_good(Sampleneuron_set(xx));
    
    nsubplot(169,169, plotplacesety{xy}, plotplacesetx{xy});
    
    title([title_for_plot{xy}]);
%     datastd = sampleneuron.datastd;
%     datamean = sampleneuron.datamean;
    All_SDFs = vertcat(Neuronlist_good(Sampleneuron_set).(SDF_for_plot{xy}));
    First_SDFs = nanmean(All_SDFs(1:2:end,:),1);
    Second_SDFs = nanmean(All_SDFs(2:2:end,:),1);
    
    plot(First_SDFs,'-b', 'LineWidth', 1);
    plot(Second_SDFs,'-r', 'LineWidth', 1);
    
%     text(400,2, ['p=' mat2str(sampleneuron.(Pvalues_in_struct{xy}),3)]);
%     text(400,2.5, ['ind=' mat2str(sampleneuron.(indices_in_struct{xy}),3)]);
    
    xlim([0,600]);
    ylim(yaxislim);
    y_lim = get(gca,'ylim');
    plot([100,100],y_lim,'k');
    %title('Population SDFs');
    ylabel('Z-scored firing rate');
end
xy = 5;
nsubplot(169,169, plotplacesety{xy}, plotplacesetx{xy});
text(0,3,'Population SDFs');
text(0,2,sprintf('n = %d', size(All_SDFs,1)/2));
ylim([0,4]);
axis off;

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(sample_neuron_folder, ['neuronpopulation.pdf']));