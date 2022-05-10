function correlation_plot_func(indices, plotpath, plotname)
StatisticalThreshold = 0.01;

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

axislabel_for_plot = {'Sensory surprise', 'Recency', 'Violation','Novelty','Reward value', 'Info anticipation'};
indices_for_plot = {'pred_vs_unpred_fam','recency_ind','violation_ind', 'pred_nov_vs_fam', 'rewardvalueindex', 'RewInfoAnticipIndex'};

flip_nov_sign = 1;

%% A simplified version just for the figure plots.
permutexy = {[4,1],[4,2],[4,3]}; 
plot_nrow = 3;
plot_ncol= 5;
plotplacesetx = {2,3,4};
plotplacesety = {1,1,1};


Include_criterion = 'noveltyexcited';% 'noveltyinhibited', 'noveltyselective', 'All'

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

for xyw = 1:length(permutexy) 
    
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
    notnanlogic = ~isnan(xaxis_ind) & ~isnan(yaxis_ind);
    
    xaxis_label = axislabel_for_plot{permutexy{xyw}(1)};
    yaxis_label = axislabel_for_plot{permutexy{xyw}(2)};
    
    
    
    %% plot a curve shows the tendency by binning
    figure(fig1);
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
    %non-related novelty neuron
    text(0.4,-0.08*(7/5), ['n = ' mat2str(sum(notnanlogic))])
    
    %axis square
    pbaspect([3,4,1]);
    % title says the criterion
    if xyw == 1
        title(Include_criterion);
    end
    
end

% save the plot
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,plotname{2}));

end