fig1 = figure;
fig2 = figure;
%Indvariablenames = {'NotNoveltySelective', 'NovelExcited', 'NovelInhibited'};
axislabel_for_plot = {'Sensory surprise', 'Recency', 'Violation','Novelty','Reward value','Reward uncertainty', 'RPE', 'Info anticipation'};
indices_for_plot = {'pred_vs_unpred_fam','recency_ind','violation_ind', 'pred_nov_vs_fam', 'rewardvalueindex','uncertaintyindex', 'rewardcuePE', 'RewInfoAnticipIndex'};
P_value_for_plot = {'Ppred_vs_unpred_fam', 'Precency_ind','Pviolation_ind', 'Ppred_nov_vs_fam', 'Prewardvalueindex', 'Puncertaintyindex', 'PrewardcuePE', 'PRewInfoAnticipIndex'};

flip_nov_sign = 1;
% Pvalues_for_plot = {'Ppred_vs_unpred_fam','Precency_ind','Pviolation_ind'};
permutexy = {[4,1],[4,2],[4,3],[4,5],[4,6],[4,7],[2,3],[1,3],[1,2], [4,8]};

plotplacesetx = {51:85,91:125,131:165,51:85,91:125,131:165,51:85,91:125,131:165, 11:45, 11:45, 11:45};
plotplacesety = {1:49, 1:49, 1:49, 61:109, 61:109, 61:109, 121:169,121:169,121:169, 61:109, 121:169, 1:49};
% barcolor = {'b','r','g'};
% Xedges = -0.5:0.05:0.5;
% Yedges = -0.5:0.05:0.5;

Include_criterion = 'noveltyexcited';% 'noveltyinhibited', 'noveltyselective'
%shuffling_num = 10000;

%Include_neurons = find(Ppred_nov_vs_fam'<StatisticalThreshold & pred_nov_vs_fam>0);
if strcmpi(Include_criterion, 'noveltyexcited')
    Include_neurons = find(Ppred_nov_vs_fam'<StatisticalThreshold & pred_nov_vs_fam>0); %1:length(pred_vs_unpred_fam);
elseif strcmpi(Include_criterion, 'noveltyinhibited')
    Include_neurons = find(Ppred_nov_vs_fam'<StatisticalThreshold & pred_nov_vs_fam<0);
elseif strcmpi(Include_criterion, 'noveltyselective')
    Include_neurons = find(Ppred_nov_vs_fam'<StatisticalThreshold);
else
    error('Invalid criterion variable');
end

NotNoveltySelective_local = ismember(Include_neurons, NotNoveltySelective);
NovelExcited_local = ismember(Include_neurons, NovelExcited);
NovelInhibited_local = ismember(Include_neurons, NovelInhibited);

for xyw = 1:length(permutexy) %%% x-axis recency ind, y-axis violation_ind
    eval(['xaxis_ind = ' indices_for_plot{permutexy{xyw}(1)} ';']);
    eval(['yaxis_ind = ' indices_for_plot{permutexy{xyw}(2)} ';']);
    
    if flip_nov_sign
        Nov_ind = pred_nov_vs_fam;
        sign_Nov = sign(Nov_ind);
        xaxis_ind = xaxis_ind.*sign_Nov;
        yaxis_ind = yaxis_ind.*sign_Nov;
    end
    
    xaxis_ind = xaxis_ind(Include_neurons);
    yaxis_ind = yaxis_ind(Include_neurons);
    
    xaxis_label = axislabel_for_plot{permutexy{xyw}(1)};
    yaxis_label = axislabel_for_plot{permutexy{xyw}(2)};
    
    figure(fig1);
    nsubplot(169,169, plotplacesety{xyw}, plotplacesetx{xyw}); set(gca,'ticklength',4*get(gca,'ticklength'))
    
    
    %% All neurons
    lim = [-1, 1];
    nbin = 25;
    binedge = linspace(lim(1),lim(2),nbin+1);
    binmid = 0.5*binedge(1:(end-1)) + 0.5*binedge(2:end);
    
    %for xyx = 1 %1: non-novelty, 2: novelty excited, 3: novelty inhibited
    %eval(['currentindset = ' Indvariablenames{xyx} ';']);
    %histogram2(xaxis_ind(currentindset),yaxis_ind(currentindset), Xedges,  Yedges,'facecolor','flat');
    [n,xedge,yedge] = histcounts2(xaxis_ind,yaxis_ind,binedge,binedge);
    
    n = n'; % transpose so that x=columns and y=rows
    
    %n(isnan(n)) = 0;
    n = n ./ nansum(n(:));
    n = log(n);
    minimum_n = min(n(n>-inf));
    n(n==-inf)=minimum_n-0.5;
    colormap(flipud(gray));
    %figuren;
    %image(binmid,binmid,colormapify(n,[0 max(n(:))],'w',interpcolor('w','k',.5),'k','w'));
    imagesc(binmid,binmid,n);
    n_contour_levels = 6;
    contour(binmid,binmid,n,n_contour_levels,'linewidth',3);
    %end
    axis([-1 1 -1 1]); %%axis square;
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    %zlabel('neuron count');
    line([-1 1],[0 0],'color',[0.3 0.3 0.3],'LineWidth',1);
    line([0 0],[-1 1],'color',[0.3 0.3 0.3],'LineWidth',1);
    %view(0,90);
    
    %%
    % now get the point which both x and y are not nan
    notnanlogic = ~isnan(xaxis_ind) & ~isnan(yaxis_ind);
    
    
%     [slope,yint] = type_2_regression(xaxis_ind(notnanlogic), yaxis_ind(notnanlogic));
%     line([-.5 .5], yint+slope.*[-0.5,0.5],'color','r','LineWidth',1);
%     %[rho,rhop] = corr(x,y,'type','Spearman');
%     
    % All datapoint value
    [rho,p] = corr(yaxis_ind(notnanlogic) , xaxis_ind(notnanlogic), 'Type', 'Spearman');
    text(0.2,0.45, ['p all ' mat2str(p,4)])
    text(0.2,0.5, ['rho all ' mat2str(rho,4)])
    % positive novelty neuron
    try
        [rho_npos,p_npos] = corr(yaxis_ind(NovelExcited_local & notnanlogic) , xaxis_ind(NovelExcited_local & notnanlogic), 'Type', 'Spearman');
        text(0.2,0.35, ['p pos ' mat2str(p_npos,4)])
        text(0.2,0.4, ['rho pos ' mat2str(rho_npos,4)])
    catch
    end
    % negative novelty neuron
    try
    [rho_nneg,p_nneg] = corr(yaxis_ind(NovelInhibited_local & notnanlogic) , xaxis_ind(NovelInhibited_local & notnanlogic), 'Type', 'Spearman');
    text(0.2,0.25, ['p neg ' mat2str(p_nneg,4)])
    text(0.2,0.3, ['rho neg ' mat2str(rho_nneg,4)])
    catch
    end
%     %non-related novelty neuron
%     [rho_nirl,p_nirl] = corr(yaxis_ind(NotNoveltySelective_local) , xaxis_ind(NotNoveltySelective_local), 'Type', 'Spearman');
    text(0.2,0.2, ['n = ' mat2str(sum(notnanlogic))])
    axis square
    % title says the criterion
    if xyw == 1
        title(Include_criterion);
    end
    
    %% plot a curve shows the tendency by binning
    figure(fig2);
    nsubplot(169,169, plotplacesety{xyw}, plotplacesetx{xyw}); set(gca,'ticklength',4*get(gca,'ticklength'))
    bins = linspace(0,1,4);%-0.5:0.2:0.5;
    %bins = linspace(-0.5,0.5,8);
    bins(1) = bins(1)-0.001;% to include the left edge
    bins(end) = bins(end)+0.001;% to include the right edge
    xdata = zeros([1,length(bins)-1]);
    ydata = zeros([1,length(bins)-1]);
    ydataerror = zeros([1,length(bins)-1]);
    for ii = 1:(length(bins)-1)
        binlogic = (xaxis_ind>=bins(ii) &  xaxis_ind<bins(ii+1));
        xdata(ii) = nanmean(xaxis_ind(binlogic));
        ydata(ii) = nanmean(yaxis_ind(binlogic));
        ydataerror(ii) = nanstd(yaxis_ind(binlogic))/sqrt(sum(~isnan(yaxis_ind(binlogic))));
        
        if sum(~isnan(yaxis_ind(binlogic)))>0
            p = signrank(yaxis_ind(binlogic));
        else
            p = nan;
        end
        
        if p<0.001
            text((bins(ii)+bins(ii+1))/2, 0.05, '***', 'FontSize', 15, 'color', 'r');
        elseif p<0.01
            text((bins(ii)+bins(ii+1))/2, 0.05, '**', 'FontSize', 15, 'color', 'r');
        elseif p<0.05
            text((bins(ii)+bins(ii+1))/2, 0.05, '*', 'FontSize', 15, 'color', 'r');
        end
    end
    
    errorbar(xdata,ydata,ydataerror,'Color','r');
    line([-1 1],[0 0],'color',[0.3 0.3 0.3],'LineWidth',1);
    line([0 0],[-1 1],'color',[0.3 0.3 0.3],'LineWidth',1);
    
    % linear fitting
    b = [ones(sum(notnanlogic),1),xaxis_ind(notnanlogic)]\[yaxis_ind(notnanlogic)];
    %objfit = fit(xaxis_ind(notnanlogic), yaxis_ind(notnanlogic), 'poly1');
    %b = [objfit.p2, objfit.p1];
    %     if strcmpi(Include_criterion, 'noveltyexcited')
    %         x_plot=0:0.1:0.5;
    %     elseif strcmpi(Include_criterion, 'noveltyinhibited')
    %         x_plot=-0.5:0.1:0;
    %     elseif strcmpi(Include_criterion, 'noveltyselective')
    %         x_plot=-0.5:0.1:0.5;
    %     else
    %         x_plot=0:0.1:0.5;
    %     end
    x_plot=0:0.1:1;
    plot(x_plot, b(1)+b(2)*x_plot, 'r');
    
    %     text(0.2,0.15, ['p remain ' mat2str(p_nirl,4)])
    %     text(0.2,0.2, ['rho remain ' mat2str(rho_nirl,4)])
    %     if strcmpi(indices_for_plot{permutexy{xyw}(1)}, 'pred_nov_vs_fam')
    %         xlim([0, 0.5]);
    %     end
    %     if strcmpi(indices_for_plot{permutexy{xyw}(2)}, 'pred_nov_vs_fam')
    %         ylim([0, 0.5]);
    %     end
    if strcmpi(Include_criterion, 'noveltyexcited')
        ylim([-0.08,0.08]*2);
    elseif strcmpi(Include_criterion, 'noveltyinhibited')
        ylim([-0.3,0.3]*2);
    elseif strcmpi(Include_criterion, 'noveltyselective')
        ylim([-0.3,0.3]*2);
    else
        ylim([-0.08,0.08]);
    end
    %ylim([-0.08,0.08]);
    xlim([-1,1]);
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    
    %% text the correlations and p values
    % All datapoint value
    [rho,p] = corr(yaxis_ind(notnanlogic) , xaxis_ind(notnanlogic), 'Type', 'Spearman');
    text(0.4,-0.08*(2/5), ['p all ' mat2str(p,4)])
    text(0.4,-0.08*(1/5), ['rho all ' mat2str(rho,4)])
    % positive novelty neuron
    try
        [rho_npos,p_npos] = corr(yaxis_ind(NovelExcited_local & notnanlogic) , xaxis_ind(NovelExcited_local & notnanlogic), 'Type', 'Spearman');
        text(0.4,-0.08*(4/5), ['p pos ' mat2str(p_npos,4)])
        text(0.4,-0.08*(3/5), ['rho pos ' mat2str(rho_npos,4)])
    catch
    end
    % negative novelty neuron
    try
    [rho_nneg,p_nneg] = corr(yaxis_ind(NovelInhibited_local & notnanlogic) , xaxis_ind(NovelInhibited_local & notnanlogic), 'Type', 'Spearman');
    text(0.4,-0.08*(6/5), ['p neg ' mat2str(p_nneg,4)])
    text(0.4,-0.08*(5/5), ['rho neg ' mat2str(rho_nneg,4)])
    catch
    end
%     %non-related novelty neuron
%     [rho_nirl,p_nirl] = corr(yaxis_ind(NotNoveltySelective_local) , xaxis_ind(NotNoveltySelective_local), 'Type', 'Spearman');
    text(0.4,-0.08*(7/5), ['n = ' mat2str(sum(notnanlogic))])
    
    axis square
    % title says the criterion
    if xyw == 1
        title(Include_criterion);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot corr(sensory surprise & recency) vs novelty index
figure(fig1);
xyw = 11;
nsubplot(169,169, plotplacesety{xyw}, plotplacesetx{xyw}); set(gca,'ticklength',4*get(gca,'ticklength'))

xaxis_ind = pred_nov_vs_fam(Include_neurons);
yaxis_ind1 = recency_ind(Include_neurons);
yaxis_ind2 = pred_vs_unpred_fam(Include_neurons);
notnanlogic = ~isnan(xaxis_ind) & ~isnan(yaxis_ind1) & ~isnan(yaxis_ind2);
xaxis_ind = xaxis_ind(notnanlogic);
yaxis_ind1 = yaxis_ind1(notnanlogic);
yaxis_ind2 = yaxis_ind2(notnanlogic);

bins = -1:0.2:1;
bins(end) = bins(end)+0.001;% to include the right edge
xdata = zeros([1,length(bins)-1]);
ydata = zeros([1,length(bins)-1]);
ylowerbound = zeros([1,length(bins)-1]);
yupperbound = zeros([1,length(bins)-1]);
for ii = 1:(length(bins)-1)
    binlogic = find(xaxis_ind>=bins(ii) &  xaxis_ind<bins(ii+1));
    if numel(binlogic)<2
        ylowerbound = nan;
        yupperbound = nan;
        xdata = nan;
        ydata = nan;
        continue;
    end
    xdata(ii) = nanmean(xaxis_ind(binlogic));
    [ydata(ii),p] = corr(yaxis_ind1(binlogic) , yaxis_ind2(binlogic), 'Type', 'Spearman');
    % ydata error: bootstrapping
    bootstrapping_corr = zeros(shuffling_num,1);
    for iii = 1:shuffling_num
        logic_sampling = randi(numel(binlogic),[numel(binlogic),1]);
        [bootstrapping_corr(iii),p] = corr(yaxis_ind1(logic_sampling) , yaxis_ind2(logic_sampling), 'Type', 'Spearman');
    end
    bootstrapping_corr = sort(bootstrapping_corr,1);
    ylowerbound(ii) = bootstrapping_corr(round(shuffling_num*0.025));
    yupperbound(ii) = bootstrapping_corr(round(shuffling_num*0.975));
end

line([-1 1],[0 0],'color',[0.3 0.3 0.3],'LineWidth',1);
errorbar(xdata,ydata,ylowerbound-ydata, yupperbound-ydata ,'Color','r');

xlabel('Novelty index')
ylabel('Correlation of recency and surprise');

text(0.2,0.2, ['n = ' mat2str(sum(notnanlogic))])
axis([-1 1 -1 1]); %%axis square;
axis square

% compare whether two corr is significantly different
figure(fig2);
xyw=12;
nsubplot(169,169, plotplacesety{xyw}, plotplacesetx{xyw}); set(gca,'ticklength',4*get(gca,'ticklength'))

Corrset1 = {[4,1],[4,2]};
Corrset2 = {[4,5],[4,6],[4,7],[4,8]};
xtext = 0;
ytext = 0;
for ii=1:numel(Corrset1)
    for jj = 1:numel(Corrset2)
        
        eval(['xaxis_ind1 = ' indices_for_plot{Corrset1{ii}(1)} '(Include_neurons);']);
        eval(['yaxis_ind1 = ' indices_for_plot{Corrset1{ii}(2)} '(Include_neurons);']);
        eval(['xaxis_ind2 = ' indices_for_plot{Corrset2{jj}(1)} '(Include_neurons);']);
        eval(['yaxis_ind2 = ' indices_for_plot{Corrset2{jj}(2)} '(Include_neurons);']);
        if flip_nov_sign
            Nov_ind = pred_nov_vs_fam(Include_neurons);
            sign_Nov = sign(Nov_ind);
            xaxis_ind1 = xaxis_ind1.*sign_Nov;
            yaxis_ind1 = yaxis_ind1.*sign_Nov;
            xaxis_ind2 = xaxis_ind2.*sign_Nov;
            yaxis_ind2 = yaxis_ind2.*sign_Nov;
        end
        
        % exclude the nan index
        excludeind = isnan(xaxis_ind1) | isnan(yaxis_ind1) | isnan(xaxis_ind2) | isnan(yaxis_ind2);
        xaxis_ind1(excludeind) = [];
        yaxis_ind1(excludeind) = [];
        xaxis_ind2(excludeind) = [];
        yaxis_ind2(excludeind) = [];
        p = Test_two_corr_difference_fun(xaxis_ind1', yaxis_ind1', xaxis_ind2', yaxis_ind2', shuffling_num);
        
        text(xtext, ytext, sprintf([axislabel_for_plot{Corrset1{ii}(1)} '+' axislabel_for_plot{Corrset1{ii}(2)} ' vs ' axislabel_for_plot{Corrset2{jj}(2)} 'p = %.5f']...
            , p), 'fontsize', 7);
        ytext = ytext+1;
    end
end
ylim([-1 ytext]);



% save the plot
figure(fig1);
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,['Indices_heatmap_all_session' '.pdf']));

% save the plot
figure(fig2);
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,['Indices_binned_errorbar_all_session' '.pdf']));

% draw 2 by 2 chart showing the number of recency-nonrecency neuron and
% surprise-nonsurprise neuron

num_nonrecency_nonsurprise = sum(Ppred_vs_unpred_fam>StatisticalThreshold & Precency_ind>StatisticalThreshold);
num_recency_nonsurprise = sum(Ppred_vs_unpred_fam>StatisticalThreshold & Precency_ind<StatisticalThreshold);
num_nonrecency_surprise = sum(Ppred_vs_unpred_fam<StatisticalThreshold & Precency_ind>StatisticalThreshold);
num_recency_surprise = sum(Ppred_vs_unpred_fam<StatisticalThreshold & Precency_ind<StatisticalThreshold);

Name = {'Non-surprise'; 'Surprise'};
Nonrecency = [num_nonrecency_nonsurprise;  num_nonrecency_surprise];
Recency = [num_recency_nonsurprise;  num_recency_surprise];
Countmatrix = [Nonrecency,Recency];
    
%T = table(Name, Nonrecency, Recency);

figure;
xyw = 1;
%nsubplot(169,169, plotplacesety{xyw}, plotplacesetx{xyw}); set(gca,'ticklength',4*get(gca,'ticklength'))
hold off
heatmap({'Non-recency', 'Recency'}, {'Non-surprise', 'Surprise'}, Countmatrix);
oddratio = Countmatrix(1,1)*Countmatrix(2,2)/(Countmatrix(1,2)*Countmatrix(2,1));
[h,p,stats] = fishertest(Countmatrix);
title(sprintf('oddratio = %0.3f, p = %f', oddratio, p));
print(gcf,'-dpdf', '-painters',fullfile(plotpath,['recency_surprise_oddratio_test ' '.pdf']));

