%barplotV2
% bar plot with indice's sign flipped whose novelty index <0

%% put recency/ surprise in a plot
%unexpected index
%unp_vs_pred_nov_odd = [Neuronlist_odd(:).unp_vs_pred_nov]'-0.5;
% NovelExcited_even = find([Neuronlist_even(:).pred_nov_vs_fam]>0.5 & [Neuronlist_even(:).P_pred_nov_vs_fam]<=StatisticalThreshold);
% NovelInhibited_even = find([Neuronlist_even(:).pred_nov_vs_fam]<0.5 & [Neuronlist_even(:).P_pred_nov_vs_fam]<=StatisticalThreshold);
% NotNoveltySelective_even = find([Neuronlist_even(:).P_pred_nov_vs_fam]>StatisticalThreshold);
%

flipsign_set = [true, false];
variables = {'pred_vs_unpred_fam', 'violation_ind','recency_ind','uncertaintyindex', 'rewardvalueindex', 'rewardcuePE', 'RewInfoAnticipIndex'};
variablenames = {'Expected vs not expected fam', 'Violation index', 'Recency index', 'Uncertainty index', 'Reward value Index', 'Reward Prediction error', 'Info anticipation index'};
Indvariablenames = {'NovelExcited', 'NovelInhibited','NotNoveltySelective', 'NoveltySelective'};
plot_positiony = {21:45, 21:45, 21:45, 21:45, 21:45, 21:45};
plot_positionx = {1:30, 41:70, 81:110, 121:150, 161:190, 195:205};
StatToCompare=0;

figure;
for xxx = 1:numel(flipsign_set)
    flipsign = flipsign_set(xxx);
    
    for xy = 1:length(Indvariablenames)
        
        nsubplot(169,209, (xxx-1)*60+plot_positiony{xy}, plot_positionx{xy}); set(gca,'ticklength',4*get(gca,'ticklength'))
        %%
        eval(['Selectcrit = ' Indvariablenames{xy} ';']);
        
        if flipsign
            datasign = sign(pred_nov_vs_fam);
        else
            datasign = ones(size(pred_nov_vs_fam));
        end
        datax = datasign(Selectcrit).*pred_vs_unpred_fam(Selectcrit);
        datay = datasign(Selectcrit).*recency_ind(Selectcrit);
        
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
            
            text(1.5,0.1, mat2str (signrank(datax, StatToCompare), 5));
            text(3.5,0.1, mat2str (signrank( datay, StatToCompare ), 5));
            text(1,-0.05, sprintf('n = %d', numel(datax)) );
            
        end
        
        ylabel('Discrimination (AUC)');
        title(Indvariablenames{xy});
        set(gca, 'xtick', [1,3], 'xticklabel', {'Sensory Surprise', 'Recency'});
        
    end
    
    
    
    p = ranksum(pred_vs_unpred_fam(NoveltySelective).*sign(pred_nov_vs_fam(NoveltySelective)), pred_vs_unpred_fam(NotNoveltySelective).*sign(pred_nov_vs_fam(NotNoveltySelective)));
    text(1,0.05, ['p comp=' mat2str(p, 4)])
    p = ranksum(recency_ind(NoveltySelective).*sign(pred_nov_vs_fam(NoveltySelective)), recency_ind(NotNoveltySelective).*sign(pred_nov_vs_fam(NotNoveltySelective)));
    text(3,0.05, ['p comp=' mat2str(p, 4)])
    %% violation index for supplemental
    %% Combine excited and inhibited neurons
    nsubplot(169,209, (xxx-1)*60+plot_positiony{5}, plot_positionx{5}); set(gca,'ticklength',4*get(gca,'ticklength'))
    %%
    if flipsign
        datasign = sign(pred_nov_vs_fam);
    else
        datasign = ones(size(pred_nov_vs_fam));
    end
    datax = datasign(NoveltySelective).*violation_ind(NoveltySelective);
    datay = datasign(NotNoveltySelective).*violation_ind(NotNoveltySelective);
    
    % datax = [violation_ind(NovelExcited); -violation_ind(NovelInhibited)];
    % datay = [violation_ind(NotNoveltySelective)];
    
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
        
        text(1.5,0.05, mat2str (signrank( datax, StatToCompare ), 5 ) )
        text(3.5,0.05, mat2str (signrank( datay, StatToCompare ), 5 ) )
        text(2.5,-0.04, mat2str (ranksum(datax, datay),4));
        %text(1,-0.05, sprintf('n = %d', numel(datax)) );
        %
        %
        %     text(1,-0.1,mat2str (( round ( ranksum( datax, dataz ) * 10000 ) ) ./ 10000 ) );
        %     text(4,-0.1,mat2str (( round ( ranksum( datay, dataz ) * 10000 ) ) ./ 10000 ) );
    end
    
    ylabel('Discrimination (AUC)');
    title('Violation');
    set(gca, 'xtick', [1,3], 'xticklabel', {'Violation/nov selective', 'Violation/other'});
    
    %text
    nsubplot(169,209, (xxx-1)*60+plot_positiony{6}, plot_positionx{6})
    p = ranksum(pred_vs_unpred_fam(NovelInhibited), pred_vs_unpred_fam(NotNoveltySelective));
    text(0,1, ['p inhibited vs other (surprise) = ' mat2str(p, 4)]);
    p = ranksum(recency_ind(NovelInhibited), recency_ind(NotNoveltySelective));
    text(0,2, ['p inhibited vs other (recency) = ' mat2str(p, 4)]);
    p = ranksum(pred_vs_unpred_fam(NovelExcited), pred_vs_unpred_fam(NotNoveltySelective));
    text(0,3, ['p excited vs other (surprise) = ' mat2str(p, 4)]);
    p = ranksum(recency_ind(NovelExcited), recency_ind(NotNoveltySelective));
    text(0,4, ['p excited vs other (recency) = ' mat2str(p, 4)]);
    
    ylim([0,5]);
    text(0,0, ['sign flipped: ' mat2str(flipsign)]);
    axis off
    
end

% save the plot
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,['Indices_barplot_all_session_V2_noveltycontrol' '.pdf']));