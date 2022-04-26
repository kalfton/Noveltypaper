function heatmapplot_func(indices, plotpath, plotname, shuffling_num)
StatisticalThreshold = 0.01;
if ~exist('shuffling_num', 'var')
    shuffling_num = 10000;
end

default_plotname = {'Indices_heatmap_all_session.pdf';
        'Indices_binned_errorbar_all_session.pdf';};
if numel(plotname)<2
    plotname(numel(plotname)+1:2) = default_plotname(numel(plotname)+1:2);
end

NovelExcited=find([indices.pred_nov_vs_fam]>0 & [indices.Ppred_nov_vs_fam]<=StatisticalThreshold)';
NovelInhibited=find([indices.pred_nov_vs_fam]<0 & [indices.Ppred_nov_vs_fam]<=StatisticalThreshold)';
NotNoveltySelective=find([[indices.Ppred_nov_vs_fam]>=StatisticalThreshold])';
NoveltySelective=find([[indices.Ppred_nov_vs_fam]<StatisticalThreshold])';

fig1 = figure;
fig2 = figure;
%Indvariablenames = {'NotNoveltySelective', 'NovelExcited', 'NovelInhibited'};
axislabel_for_plot = {'Sensory surprise', 'Recency', 'Violation','Novelty','Reward value', 'Info anticipation'};
indices_for_plot = {'pred_vs_unpred_fam','recency_ind','violation_ind', 'pred_nov_vs_fam', 'rewardvalueindex', 'RewInfoAnticipIndex'};
%P_value_for_plot = {'Ppred_vs_unpred_fam', 'Precency_ind','Pviolation_ind', 'Ppred_nov_vs_fam', 'Prewardvalueindex', 'Puncertaintyindex', 'PrewardcuePE', 'PRewInfoAnticipIndex'};

flip_nov_sign = 1;
% Pvalues_for_plot = {'Ppred_vs_unpred_fam','Precency_ind','Pviolation_ind'};

%% A simplified version just for the figure plots.
permutexy = {[4,1],[4,2],[4,3], [1,2]}; 
plot_nrow = 3;
plot_ncol= 5;
plotplacesetx = {2,3,4,2,1};
plotplacesety = {1,1,1,2,1};


Include_criterion = 'noveltyexcited';% 'noveltyinhibited', 'noveltyselective', 'All'

%Include_neurons = find(Ppred_nov_vs_fam'<StatisticalThreshold & indices.pred_nov_vs_fam>0);
if strcmpi(Include_criterion, 'noveltyexcited')
    Include_neurons = find(indices.Ppred_nov_vs_fam<StatisticalThreshold & indices.pred_nov_vs_fam>0); %1:length(pred_vs_unpred_fam);
elseif strcmpi(Include_criterion, 'noveltyinhibited')
    Include_neurons = find(indices.Ppred_nov_vs_fam<StatisticalThreshold & indices.pred_nov_vs_fam<0);
elseif strcmpi(Include_criterion, 'noveltyselective')
    Include_neurons = find(indices.Ppred_nov_vs_fam<StatisticalThreshold);
elseif strcmpi(Include_criterion, 'All')
    Include_neurons = find(indices.Ppred_nov_vs_fam<inf);
    flip_nov_sign = 0;
else
    error('Invalid criterion variable');
end

NotNoveltySelective_local = ismember(Include_neurons, NotNoveltySelective);
NovelExcited_local = ismember(Include_neurons, NovelExcited);
NovelInhibited_local = ismember(Include_neurons, NovelInhibited);

for xyw = 1:length(permutexy) %%% x-axis recency ind, y-axis violation_ind
    
    xaxis_ind = indices.(indices_for_plot{permutexy{xyw}(1)});
    yaxis_ind = indices.(indices_for_plot{permutexy{xyw}(2)});
    
    if flip_nov_sign
        Nov_ind = indices.pred_nov_vs_fam;
        sign_Nov = sign(Nov_ind);
        xaxis_ind = xaxis_ind.*sign_Nov;
        yaxis_ind = yaxis_ind.*sign_Nov;
    end
    
    xaxis_ind = xaxis_ind(Include_neurons);
    yaxis_ind = yaxis_ind(Include_neurons);
    
    xaxis_label = axislabel_for_plot{permutexy{xyw}(1)};
    yaxis_label = axislabel_for_plot{permutexy{xyw}(2)};
    
    figure(fig1);
    nsubplot(plot_nrow, plot_ncol, plotplacesety{xyw}, plotplacesetx{xyw}); set(gca,'ticklength',4*get(gca,'ticklength'));
    
    
    %% All neurons
    lim = [-1, 1];
    nbin = 25;
    binedge = linspace(lim(1),lim(2),nbin+1);
    binmid = 0.5*binedge(1:(end-1)) + 0.5*binedge(2:end);
    
    [n,~,~] = histcounts2(xaxis_ind,yaxis_ind,binedge,binedge);
    
    n = n'; % transpose so that x=columns and y=rows
    
    %n(isnan(n)) = 0;
    n = n ./ nansum(n(:));
    n = log(n);
    minimum_n = min(n(n>-inf));
    n(n==-inf)=minimum_n-0.5;
    colormap(flipud(gray));
    %figuren;
    imagesc(binmid,binmid,n);
    n_contour_levels = 6;
    contour(binmid,binmid,n,n_contour_levels,'linewidth',3);
    %end
    axis([-1 1 -1 1]); %%axis square;
    xlabel(xaxis_label)
    ylabel(yaxis_label)
    line([-1 1],[0 0],'color',[0.3 0.3 0.3],'LineWidth',1);
    line([0 0],[-1 1],'color',[0.3 0.3 0.3],'LineWidth',1);
    
    %%
    % now get the point which both x and y are not nan
    notnanlogic = ~isnan(xaxis_ind) & ~isnan(yaxis_ind);
    
    
%     
    % All datapoint value
    [rho,p] = corr(yaxis_ind(notnanlogic) , xaxis_ind(notnanlogic), 'Type', 'Spearman');
    text(0.2,0.7, ['p all ' mat2str(p,4)])
    text(0.2,0.8, ['rho all ' mat2str(rho,4)])
    % positive novelty neuron
    try
        [rho_npos,p_npos] = corr(yaxis_ind(NovelExcited_local & notnanlogic) , xaxis_ind(NovelExcited_local & notnanlogic), 'Type', 'Spearman');
        text(0.2,0.5, ['p pos ' mat2str(p_npos,4)])
        text(0.2,0.6, ['rho pos ' mat2str(rho_npos,4)])
    catch
    end
    % negative novelty neuron
    try
    [rho_nneg,p_nneg] = corr(yaxis_ind(NovelInhibited_local & notnanlogic) , xaxis_ind(NovelInhibited_local & notnanlogic), 'Type', 'Spearman');
    text(0.2,0.3, ['p neg ' mat2str(p_nneg,4)])
    text(0.2,0.4, ['rho neg ' mat2str(rho_nneg,4)])
    catch
    end
%     %non-related novelty neuron
    text(0.2,0.2, ['n = ' mat2str(sum(notnanlogic))])
    axis square
    % title says the criterion
    if xyw == 1
        title(Include_criterion);
    end
    
    %% plot a curve shows the tendency by binning
    figure(fig2);
    nsubplot(plot_nrow, plot_ncol, plotplacesety{xyw}, plotplacesetx{xyw}); set(gca,'ticklength',4*get(gca,'ticklength'));
    bins = linspace(0,1,4);%-0.5:0.2:0.5;
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
    x_plot=0:0.1:1;
    plot(x_plot, b(1)+b(2)*x_plot, 'r');
    
    if strcmpi(Include_criterion, 'noveltyexcited')
        ylim([-0.08,0.08]*2);
    elseif strcmpi(Include_criterion, 'noveltyinhibited')
        ylim([-0.3,0.3]*2);
    elseif strcmpi(Include_criterion, 'noveltyselective')
        ylim([-0.3,0.3]*2);
    else
        ylim([-0.08,0.08]);
    end
    if contains(indices_for_plot{permutexy{xyw}(2)},'violation')
        ylim([-0.04,0.16])
        set(gca, 'ytick', -0.04:0.04:0.16);
    else
        ylim([-0.04,0.12])
        set(gca, 'ytick', -0.04:0.04:0.16);
    end
    %ylim([-0.08,0.08]);
    xlim([0,1]);
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
    text(0.4,-0.08*(7/5), ['n = ' mat2str(sum(notnanlogic))])
    
    %axis square
    pbaspect([3,4,1]);
    % title says the criterion
    if xyw == 1
        title(Include_criterion);
    end
    
end

%% plot corr(sensory surprise & recency) vs novelty index
if ~strcmpi(Include_criterion, 'All')

% compare whether two corr is significantly different
figure(fig2);
nsubplot(plot_nrow, plot_ncol, plotplacesety{end}, plotplacesetx{end}); set(gca,'ticklength',4*get(gca,'ticklength'))

Corrset1 = {[4,1],[4,2]};
Corrset2 = {[4,5],[4,6]};
xtext = 0;
ytext = 0;
for ii=1:numel(Corrset1)
    for jj = 1:numel(Corrset2)
        xaxis_ind1 = indices.(indices_for_plot{Corrset1{ii}(1)})(Include_neurons);
        yaxis_ind1 = indices.(indices_for_plot{Corrset1{ii}(2)})(Include_neurons);
        xaxis_ind2 = indices.(indices_for_plot{Corrset2{jj}(1)})(Include_neurons);
        yaxis_ind2 = indices.(indices_for_plot{Corrset2{jj}(2)})(Include_neurons);
        
        if flip_nov_sign
            Nov_ind = indices.pred_nov_vs_fam(Include_neurons);
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
end


% save the plot
figure(fig1);
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,plotname{1}));

% save the plot
figure(fig2);
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,plotname{2}));

end