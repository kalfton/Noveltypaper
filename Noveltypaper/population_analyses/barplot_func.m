function barplot_func(indices, plotpath, plotname)
StatisticalThreshold = 0.01;
NovelExcited=find([indices.pred_nov_vs_fam]>0 & [indices.Ppred_nov_vs_fam]<=StatisticalThreshold)';
NovelInhibited=find([indices.pred_nov_vs_fam]<0 & [indices.Ppred_nov_vs_fam]<=StatisticalThreshold)';
NotNoveltySelective=find([[indices.Ppred_nov_vs_fam]>=StatisticalThreshold])';
NoveltySelective=find([[indices.Ppred_nov_vs_fam]<StatisticalThreshold])';

flipsign_set = [true, false];
Indvariablenames = {'NovelExcited', 'NovelInhibited','NotNoveltySelective', 'NoveltySelective'};
plot_nrow = 3;
plot_ncol= 6;
plot_positiony = {1, 1, 1, 1, 1, 1};
plot_positionx = {1, 2, 3, 4, 5, 6};
StatToCompare=0;

figure;
for xxx = 1:numel(flipsign_set)
    flipsign = flipsign_set(xxx);
    
    for xy = 1:length(Indvariablenames)
        
        nsubplot(plot_nrow, plot_ncol, (xxx-1)+plot_positiony{xy}, plot_positionx{xy}); set(gca,'ticklength',4*get(gca,'ticklength'))
        %%
        eval(['Selectcrit = ' Indvariablenames{xy} ';']);
        
        if flipsign
            datasign = sign(indices.pred_nov_vs_fam);
        else
            datasign = ones(size(indices.pred_nov_vs_fam));
        end
        datax = datasign(Selectcrit).*indices.pred_vs_unpred_fam(Selectcrit);
        datay = datasign(Selectcrit).*indices.recency_ind(Selectcrit);
        
        datax = datax(~isnan(datax));
        datay = datay(~isnan(datay));
        
        %%
        if ~isempty(datax) && ~isempty(datay)
            bar(1,mean(datax,1),'w')
            bar(3,mean(datay,1),'w')
            ylim([-.1  .1])
            %
            JitterVar=0.2; %jitter for scatter plot
            R=datax;
            R(:,2)=1.5;
            %scatter(R(:,2),R(:,1),10,'r','filled','jitter','on', 'jitterAmount',JitterVar); hold on;
            %
            R=datay;
            R(:,2)=3.5;
            %scatter(R(:,2),R(:,1),10,'b','filled','jitter','on', 'jitterAmount',JitterVar); hold on;
            %
            errorbar([1],mean(datax,1), std(datax,1)./sqrt(size(datax,1)),'k','LineWidth',2)
            errorbar([3],mean(datay,1), std(datay,1)./sqrt(size(datay,1)),'k','LineWidth',2)
            %
            
            text(1.5,0.1, sprintf('p= %.5f', signrank(datax, StatToCompare)));
            text(3.5,0.1, sprintf('p= %.5f', signrank( datay, StatToCompare )));
            text(1,-0.05, sprintf('n = %d', numel(datax)) );
            
            text(1.5,0.15, sprintf('ave= %.5f', mean(datax)));
            text(3.5,0.15, sprintf('ave= %.5f', mean(datay)));
        end
        
        ylabel('Discrimination (AUC)');
        title(Indvariablenames{xy});
        set(gca, 'xtick', [1,3], 'xticklabel', {'Sensory Surprise', 'Recency'});
        
    end
    
    
    
    p = ranksum(indices.pred_vs_unpred_fam(NoveltySelective).*sign(indices.pred_nov_vs_fam(NoveltySelective)), indices.pred_vs_unpred_fam(NotNoveltySelective).*sign(indices.pred_nov_vs_fam(NotNoveltySelective)));
    text(1,0.05, ['p comp=' mat2str(p, 4)])
    p = ranksum(indices.recency_ind(NoveltySelective).*sign(indices.pred_nov_vs_fam(NoveltySelective)), indices.recency_ind(NotNoveltySelective).*sign(indices.pred_nov_vs_fam(NotNoveltySelective)));
    text(3,0.05, ['p comp=' mat2str(p, 4)])
    %% violation index for supplemental
    %% Combine excited and inhibited neurons
    nsubplot(plot_nrow, plot_ncol, (xxx-1)+plot_positiony{5}, plot_positionx{5}); set(gca,'ticklength',4*get(gca,'ticklength'))
    %%
    if flipsign
        datasign = sign(indices.pred_nov_vs_fam);
    else
        datasign = ones(size(indices.pred_nov_vs_fam));
    end
    datax = datasign(NoveltySelective).*indices.violation_ind(NoveltySelective);
    datay = datasign(NotNoveltySelective).*indices.violation_ind(NotNoveltySelective);
    
    % datax = [indices.violation_ind(NovelExcited); -indices.violation_ind(NovelInhibited)];
    % datay = [indices.violation_ind(NotNoveltySelective)];
    
    datax = datax(~isnan(datax));
    datay = datay(~isnan(datay));
    
    %%
    if ~isempty(datax) && ~isempty(datay)
        bar(1,mean(datax,1),'w')
        bar(3,mean(datay,1),'w')
        ylim([-.1  .1])
        %
        errorbar([1],mean(datax,1), std(datax,1)./sqrt(size(datax,1)),'k','LineWidth',2)
        errorbar([3],mean(datay,1), std(datay,1)./sqrt(size(datay,1)),'k','LineWidth',2)
        %
        
        text(1.5,0.05, mat2str (signrank( datax, StatToCompare ), 5 ) )
        text(3.5,0.05, mat2str (signrank( datay, StatToCompare ), 5 ) )
        text(2.5,-0.04, mat2str (ranksum(datax, datay),4));
    end
    
    ylabel('Discrimination (AUC)');
    title('Violation');
    set(gca, 'xtick', [1,3], 'xticklabel', {'Violation/nov selective', 'Violation/other'});
    
    %text
    nsubplot(plot_nrow, plot_ncol, (xxx-1)+plot_positiony{6}, plot_positionx{6});
    p = ranksum(indices.pred_vs_unpred_fam(NovelInhibited), indices.pred_vs_unpred_fam(NotNoveltySelective));
    text(0,0, ['p inhibited vs other (surprise) = ' mat2str(p, 4)]);
    p = ranksum(indices.recency_ind(NovelInhibited), indices.recency_ind(NotNoveltySelective));
    text(0,1, ['p inhibited vs other (recency) = ' mat2str(p, 4)]);
    p = ranksum(indices.pred_vs_unpred_fam(NovelExcited), indices.pred_vs_unpred_fam(NotNoveltySelective));
    text(0,2, ['p excited vs other (surprise) = ' mat2str(p, 4)]);
    p = ranksum(indices.recency_ind(NovelExcited), indices.recency_ind(NotNoveltySelective));
    text(0,3, ['p excited vs other (recency) = ' mat2str(p, 4)]);
    
    ylim([0,5]);
    text(0,4, ['sign flipped: ' mat2str(flipsign)], 'FontSize', 14, 'Color', 'r');
    axis off
    
end

% save the plot
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,plotname));

end