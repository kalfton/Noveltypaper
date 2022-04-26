%%% bar plots
figure;
variables = {'pred_vs_unpred_fam', 'violation_ind','recency_ind','uncertaintyindex', 'rewardvalueindex', 'rewardcuePE', 'RewInfoAnticipIndex'};
variablenames = {'Expected vs not expected fam', 'Violation index', 'Recency index', 'Uncertainty index', 'Reward value Index', 'Reward Prediction error', 'Info anticipation index'};
plot_positiony = {1:25, 41:65, 81:105, 1:25, 41:65, 81:105, 1:25};
plot_positionx = {1:30,1:30,1:30, 41:70,41:70,41:70, 81:110};
StatToCompare=0;

for xy = 1:3;%length(variables)

nsubplot(169,209, plot_positiony{xy}, plot_positionx{xy}); set(gca,'ticklength',4*get(gca,'ticklength'))
%%
eval(['IndexToPlot = ' variables{xy} ';']);

datax= [IndexToPlot(NovelExcited);];
datay= [IndexToPlot(NovelInhibited)];
dataz= IndexToPlot(NotNoveltySelective);
datax = datax(~isnan(datax));
datay = datay(~isnan(datay));
dataz = dataz(~isnan(dataz));
%%
if ~isempty(datax) && ~isempty(datay) && ~isempty(dataz)
    bar(1,mean(datax,1),'w')
    bar(3,mean(datay,1),'w')
    bar(5,mean(dataz,1),'w')
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
    R=dataz;
    R(:,2)=5.5;
    %scatter(R(:,2),R(:,1),10,'k','filled','jitter','on', 'jitterAmount',JitterVar); hold on;
    %
    errorbar([1],mean(datax,1), std(datax,1)./sqrt(size(datax,1)),'k','LineWidth',2)
    errorbar([3],mean(datay,1), std(datay,1)./sqrt(size(datay,1)),'k','LineWidth',2)
    errorbar([5],mean(dataz,1), std(dataz,1)./sqrt(size(dataz,1)),'k','LineWidth',2)
    %
    
    text(1.5,0.05, mat2str (( round ( signrank( datax, StatToCompare ) * 10000 ) ) ./ 10000 ) )
    text(3.5,0.05, mat2str (( round ( signrank( datay, StatToCompare ) * 10000 ) ) ./ 10000 ) )
    text(5.5,0.05, mat2str (( round ( signrank( dataz, StatToCompare ) * 10000 ) ) ./ 10000 ) )
    %
    %
    text(1,-0.04,mat2str (( round ( ranksum( datax, dataz ) * 10000 ) ) ./ 10000 ) );
    text(4,-0.04,mat2str (( round ( ranksum( datay, dataz ) * 10000 ) ) ./ 10000 ) );
end

ylabel('Discrimination (AUC)')
xlabel(variablenames{xy})
set(gca, 'xtick', [1,3,5], 'xticklabel', {'Nov Excited', 'Nov Inhibited', 'Other'});

end


%% put recency/ surprise in a plot

variables = {'pred_vs_unpred_fam', 'violation_ind','recency_ind','uncertaintyindex', 'rewardvalueindex', 'rewardcuePE', 'RewInfoAnticipIndex'};
variablenames = {'Expected vs not expected fam', 'Violation index', 'Recency index', 'Uncertainty index', 'Reward value Index', 'Reward Prediction error', 'Info anticipation index'};
Indvariablenames = {'NovelExcited', 'NovelInhibited','NotNoveltySelective'};
plot_positiony = {41:65, 41:65, 41:65, 41:65, 81:105, 1:25};
plot_positionx = {81:110, 121:150, 161:190, 41:70, 41:70, 41:70};
StatToCompare=0;

for xy = 1:3%length(variables)

nsubplot(169,209, plot_positiony{xy}, plot_positionx{xy}); set(gca,'ticklength',4*get(gca,'ticklength'))
%%
eval(['Selectcrit = ' Indvariablenames{xy} ';']);

datax = pred_vs_unpred_fam(Selectcrit);
datay = recency_ind(Selectcrit);

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
    
    text(1.5,0.1, mat2str (( round ( signrank( datax, StatToCompare ) * 10000 ) ) ./ 10000 ) )
    text(3.5,0.1, mat2str (( round ( signrank( datay, StatToCompare ) * 10000 ) ) ./ 10000 ) )
    text(1,-0.05, sprintf('n = %d', numel(datax)) );
    %
    %
%     text(1,-0.1,mat2str (( round ( ranksum( datax, dataz ) * 10000 ) ) ./ 10000 ) );
%     text(4,-0.1,mat2str (( round ( ranksum( datay, dataz ) * 10000 ) ) ./ 10000 ) );
end

ylabel('Discrimination (AUC)');
title(Indvariablenames{xy});
set(gca, 'xtick', [1,3], 'xticklabel', {'Sensory Surprise', 'Recency'});

end

ranksum(pred_vs_unpred_fam(NovelInhibited), pred_vs_unpred_fam(NotNoveltySelective))
ranksum(recency_ind(NovelInhibited), recency_ind(NotNoveltySelective))


%% Combine excited and inhibited neurons
nsubplot(169,209, plot_positiony{4}, plot_positionx{4}); set(gca,'ticklength',4*get(gca,'ticklength'))
%%

datax = [pred_vs_unpred_fam(NovelExcited); -pred_vs_unpred_fam(NovelInhibited)];
datay = [recency_ind(NovelExcited); -recency_ind(NovelInhibited)];

datax = datax(~isnan(datax));
datay = datay(~isnan(datay));

%%
if ~isempty(datax) && ~isempty(datay) && ~isempty(dataz)
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
    
    text(1.5,0.1, mat2str (( round ( signrank( datax, StatToCompare ) * 10000 ) ) ./ 10000 ) )
    text(3.5,0.1, mat2str (( round ( signrank( datay, StatToCompare ) * 10000 ) ) ./ 10000 ) )
    text(1,-0.05, sprintf('n = %d', numel(datax)) );
    %
    %
%     text(1,-0.1,mat2str (( round ( ranksum( datax, dataz ) * 10000 ) ) ./ 10000 ) );
%     text(4,-0.1,mat2str (( round ( ranksum( datay, dataz ) * 10000 ) ) ./ 10000 ) );
end

ylabel('Discrimination (AUC)');
title('Novelty selective');
set(gca, 'xtick', [1,3], 'xticklabel', {'Sensory Surprise', 'Recency'});

p = ranksum([pred_vs_unpred_fam(NovelExcited); -pred_vs_unpred_fam(NovelInhibited)], pred_vs_unpred_fam(NotNoveltySelective));
text(1,0.05, ['p comp=' mat2str(p, 4)])
p = ranksum([recency_ind(NovelExcited); -recency_ind(NovelInhibited)], recency_ind(NotNoveltySelective));
text(3,0.05, ['p comp=' mat2str(p, 4)])
%% violation index for supplemental
%% Combine excited and inhibited neurons
nsubplot(169,209, plot_positiony{5}, plot_positionx{5}); set(gca,'ticklength',4*get(gca,'ticklength'))
%%

datax = [violation_ind(NovelExcited); -violation_ind(NovelInhibited)];
datay = [violation_ind(NotNoveltySelective)];

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
    
    text(1.5,0.05, mat2str (( round ( signrank( datax, StatToCompare ) * 10000 ) ) ./ 10000 ) )
    text(3.5,0.05, mat2str (( round ( signrank( datay, StatToCompare ) * 10000 ) ) ./ 10000 ) )
    text(2.5,-0.04, mat2str (ranksum(datax, datay),4));
    %text(1,-0.05, sprintf('n = %d', numel(datax)) );
    %
    %
%     text(1,-0.1,mat2str (( round ( ranksum( datax, dataz ) * 10000 ) ) ./ 10000 ) );
%     text(4,-0.1,mat2str (( round ( ranksum( datay, dataz ) * 10000 ) ) ./ 10000 ) );
end

ylabel('Discrimination (AUC)');
title('Novelty selective');
set(gca, 'xtick', [1,3], 'xticklabel', {'Violation/nov selective', 'Violation/other'});


% save the plot
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',fullfile(plotpath,['Indices_barplot_all_session' '.pdf']));