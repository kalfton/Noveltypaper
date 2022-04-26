% time scale analysis

x_index = {'recency_ind_match_pos'};
y_index = {'withindaylearning_newindex', 'acrossdayforget_newindex'};
Include_neurons = [Neuronlist_good.learningforgetinganalysis];% & [Neuronlist_good.intrinsic_time]<5000 & [Neuronlist_good.intrinsic_time]>10;

plotplacesetx = {51:85,91:125,131:165,51:85,91:125,131:165,51:85,91:125,131:165, 11:45, 11:45, 11:45};
plotplacesety = {1:49, 1:49, 1:49, 61:109, 61:109, 61:109, 121:169,121:169,121:169, 61:109, 121:169, 1:49};

figure();
for ii = 1:numel(y_index) %%% x-axis recency ind, y-axis violation_ind
%     eval(['xaxis_ind = ' indices_for_plot{permutexy{xyw}(1)} ';']);
%     eval(['yaxis_ind = ' indices_for_plot{permutexy{xyw}(2)} ';']);
    
    xaxis_ind = [Neuronlist_good(Include_neurons).(x_index{1})]';
    yaxis_ind = [Neuronlist_good(Include_neurons).(y_index{ii})]';
    
    
    xaxis_label = x_index{1};
    yaxis_label = y_index{ii};
    
    nsubplot(169,169, plotplacesety{ii}, plotplacesetx{ii}); set(gca,'ticklength',4*get(gca,'ticklength'))
    
    scatter(xaxis_ind,yaxis_ind,2, 'filled');
%     %% All neurons
%     lim = [-1, 1];
%     nbin = 25;
%     binedge = linspace(lim(1),lim(2),nbin+1);
%     binmid = 0.5*binedge(1:(end-1)) + 0.5*binedge(2:end);
%     
%     [n,xedge,yedge] = histcounts2(xaxis_ind,yaxis_ind,binedge,binedge);
%     
%     n = n'; % transpose so that x=columns and y=rows
%     
%     %n(isnan(n)) = 0;
%     n = n ./ nansum(n(:));
%     n = log(n);
%     minimum_n = min(n(n>-inf));
%     n(n==-inf)=minimum_n-0.5;
%     colormap(flipud(gray));
%     %figuren;
%     %image(binmid,binmid,colormapify(n,[0 max(n(:))],'w',interpcolor('w','k',.5),'k','w'));
%     imagesc(binmid,binmid,n);
%     n_contour_levels = 6;
%     contour(binmid,binmid,n,n_contour_levels,'linewidth',3);
%     %end
%     axis([-1 1 -1 1]); %%axis square;
    xlabel(xaxis_label, 'interpreter', 'None');
    ylabel(yaxis_label, 'interpreter', 'None')
    %zlabel('neuron count');
    line([-1 1],[0 0],'color',[0.3 0.3 0.3],'LineWidth',1);
    line([0 0],[-1 1],'color',[0.3 0.3 0.3],'LineWidth',1);
    %view(0,90);
    
    %%
    % now get the point which both x and y are not nan
    notnanlogic = ~isnan(xaxis_ind) & ~isnan(yaxis_ind);
    
    % All datapoint value
    [rho,p] = corr(yaxis_ind(notnanlogic) , xaxis_ind(notnanlogic), 'Type', 'Spearman');
    text(0.2,0.7, ['p all ' mat2str(p,4)])
    text(0.2,0.8, ['rho all ' mat2str(rho,4)])
    
%     %non-related novelty neuron
%     [rho_nirl,p_nirl] = corr(yaxis_ind(NotNoveltySelective_local) , xaxis_ind(NotNoveltySelective_local), 'Type', 'Spearman');
    text(0.2,0.2, ['n = ' mat2str(sum(notnanlogic))])
    axis square
    
end

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,'timescale_correlation_all.pdf'));