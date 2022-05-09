function pupil_plots(pupillist, Monkey, plotpath)
% pupil plots
% variables = {'pred_nov_vs_fam','pred_vs_unpred_fam','recency_ind_match_pos'};
variablenames = {'Novelty index', 'Sensory surprise index', 'Recency index'};
plot_nrow = 2;
plot_ncol= 4;
plot_positiony = {1, 1, 1, 1};
plot_positionx = {1, 2, 3, 4};

StatToCompare=0;
varthr = 0.001;

includeind = zeros(size(pupillist));
for ii  = 1:numel(pupillist)
    if any(isnan(pupillist(ii).pred_nov_vs_fam_sdfs(:))) || any(isnan(pupillist(ii).pred_vs_unpred_fam_sdfs(:))) || any(isnan(pupillist(ii).recency_ind__match_pos_sdfs(:)))
        includeind(ii) = 0;
    elseif var(pupillist(ii).pred_nov_vs_fam_sdfs(:))<varthr || var(pupillist(ii).pred_vs_unpred_fam_sdfs(:))<varthr || var(pupillist(ii).recency_ind__match_pos_sdfs(:))<varthr
        includeind(ii) = 0;
    else
        includeind(ii) = 1;
    end
end
includeind = logical(includeind);

%plot
figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation

nsubplot(plot_nrow, plot_ncol, plot_positiony{1}, plot_positionx{1}); set(gca,'ticklength',4*get(gca,'ticklength'))

datax= -([pupillist.pred_nov_vs_fam]'*2-1);
datay= -([pupillist.pred_vs_unpred_fam]'*2-1);
dataz= -([pupillist.recency_ind_match_pos]'*2-1);

datax= datax(includeind);
datay= datay(includeind);
dataz= dataz(includeind);

bar(1,mean(datax,1),'w')
bar(3,mean(datay,1),'w')
bar(5,mean(dataz,1),'w')
ylim([-.4 .6])
%
errorbar([1],mean(datax,1), std(datax,1)./sqrt(size(datax,1)),'k','LineWidth',2)
errorbar([3],mean(datay,1), std(datay,1)./sqrt(size(datay,1)),'k','LineWidth',2)
errorbar([5],mean(dataz,1), std(dataz,1)./sqrt(size(dataz,1)),'k','LineWidth',2)
%

text(1.5,0.1, mat2str (( round ( signrank( datax, StatToCompare ) * 10000 ) ) ./ 10000 ) )
text(3.5,0.1, mat2str (( round ( signrank( datay, StatToCompare ) * 10000 ) ) ./ 10000 ) )
text(5.5,0.1, mat2str (( round ( signrank( dataz, StatToCompare ) * 10000 ) ) ./ 10000 ) )
%

text(4,-0.1,['n = ', mat2str(numel(datax))]);

ylabel('Pupil contraction (-pupil diameter/area)')
set(gca, 'xtick', [1,3,5], 'xticklabel', variablenames, 'xticklabelRotation', -45);



%% SDF plot of pupil data

% the PDS pupil data has dentate shape and here I use a low pass filter
% to get rid of it
[b_,a_] = butter(6,0.1*2); % cutoff frequency at 2*pi*0.1 radius per sample

sdf_variables = {'pred_nov_vs_fam_sdfs',  'pred_vs_unpred_fam_sdfs', 'recency_ind__match_pos_sdfs'};

sdf_legend = {{'familiar', 'novel'},{'expected', 'surprising'},{'recent','nonrecent'}};

for ii = 1:numel(sdf_variables)
    nsubplot(plot_nrow, plot_ncol, plot_positiony{ii+1}, plot_positionx{ii+1}); set(gca,'ticklength',4*get(gca,'ticklength'))
    pupilsdfs = vertcat(pupillist(includeind).(sdf_variables{ii}));
    pupilsdfs = filtfilt(b_, a_, pupilsdfs')'; % filt the sdfs.
    pupilsdfs_1 = pupilsdfs(1:2:end,:);
    pupilsdfs_2 = pupilsdfs(2:2:end,:);
    plot_x = 1:size(pupilsdfs,2);
    
    try
        plt1 = shadedErrorBar(plot_x,pupilsdfs_1,{@(x) nanmean(x,1) ,@(x) nanstd(x,1)./sqrt(size(x,1))},{'-r', 'LineWidth', 0.5}, 0);
        plt2 = shadedErrorBar(plot_x,pupilsdfs_2,{@(x) nanmean(x,1) ,@(x) nanstd(x,1)./sqrt(size(x,1))},{'-b', 'LineWidth', 0.5}, 0);
    catch
        plot(plot_x,mean(pupilsdfs_1,1), '-r');
        plot(plot_x,mean(pupilsdfs_2,1), '-b');
    end
    
    ylim([-1,1]);
    xlim([0,600]);
    plot([100,100], ylim, 'k');
    title(variablenames{ii});
    ylabel('Pupil diameter/area');
    
    lgd = legend([plt1.mainLine, plt2.mainLine], sdf_legend{ii}, 'interpreter', 'None');
end

print(gcf,'-dpdf', '-painters',fullfile(plotpath, ['Indices_barplot_and_sdfs_pupil_NFL_Monkey_' Monkey '.pdf']));

end
