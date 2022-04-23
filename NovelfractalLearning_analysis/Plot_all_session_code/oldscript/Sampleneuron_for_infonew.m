%%%% find a good reward and info encoding neuron and plot its SDF

Sampleneuron_set = find([Neuronlist_good(:).rewardvalueindex_precue]>0 & [Neuronlist_good(:).rewardvalueindexP_precue]<StatisticalThreshold &...
     [Neuronlist_good(:).RewInfoAnticipIndex_split]>0 & [Neuronlist_good(:).RewInfoAnticipIndexP_split]<StatisticalThreshold);
yaxislim = [-1,3];
id =6;
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

indices_in_struct = {'rewardvalueindex_precue', 'RewInfoAnticipIndex_split', 'RewInfoAnticipIndex_split'};
Pvalues_in_struct = {'rewardvalueindexP_precue', 'RewInfoAnticipIndexP_split', 'RewInfoAnticipIndexP_split'};
sdf_legend = {{'noreward info', 'reward info'},{'certain info', 'uncertain info'}, {'certain noinfo', 'uncertain noinfo'}};
title_for_plot = {'Reward value', 'Info anticipation info', 'Info anticipation noinfo'};

sdfname = {'SDF_rew100ni'
'SDF_rew50ni_del'
'SDF_rew50ni_ndel'
'SDF_rew0ni'
'SDF_rew100i'
'SDF_rew50i_del'
'SDF_rew50i_ndel'
'SDF_rew0i'};
sdf_ind_pair = { {[8,4],[5,1]}, {[8,5],[6,7]}, {[1,4],[2,3]} };

plotplacesetx = {10:30, 40:60,70:90,100:120,130:150};
plotplacesety = {5:25, 5:25, 5:25, 5:25, 5:25};
ygap = 30;

figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
figurenum  = 1;
rownum = 0;
for xx = 1:length(Sampleneuron_set)
    if (rownum+1)*ygap>169
        % save current figure
        print(gcf,'-dpdf', '-painters',fullfile(['sampleneurons_info' mat2str(figurenum) '.pdf']));
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
        
        current_pair = sdf_ind_pair{xy};
        indset1 = current_pair{1};
        indset2 = current_pair{2};
        
        %Extract the sdf of selective neurons here.
        sdfs_1 = [];
        sdfs_2 = [];
        for iii = 1:numel(indset1)
            tempsdfs = vertcat(sampleneuron.(sdfname{indset1(iii)}));
            sdfs_1 = [sdfs_1; tempsdfs];
        end
        for iii = 1:numel(indset2)
            tempsdfs = vertcat(sampleneuron.(sdfname{indset2(iii)}));
            sdfs_2 = [sdfs_2; tempsdfs];
        end
        %n_neuron = sum(sig_logi);
        
        
        plot(datamean+datastd*mean(sdfs_1,1),'-b', 'LineWidth', 1);
        plot(datamean+datastd*mean(sdfs_2,1),'-r', 'LineWidth', 1);
        if xy==1
            localyaxislim = yaxislim*datastd+datamean;
        end
        
        text(1800,2, ['p=' mat2str(sampleneuron.(Pvalues_in_struct{xy}),3)]);
        text(1800,8, ['ind=' mat2str(sampleneuron.(indices_in_struct{xy}),3)]);
        
        xlim([0,2000]);
        ylim(localyaxislim);
        y_lim = get(gca,'ylim');
        plot([1000,1000],y_lim,'k');
    end
    % electrode ID and filename
    xy = 5;
    nsubplot(169,169, rownum*ygap+plotplacesety{xy}, plotplacesetx{xy});
    text(0,1,sampleneuron.filename);
    text(0,2,['Depth ' mat2str(sampleneuron.depth,3)]);
    text(0,3,sampleneuron.name);
    text(0,4,char(sampleneuron.region));
    ylim([0,5]);
    axis off;
    
    rownum = rownum+1;
end
print(gcf,'-dpdf', '-painters',fullfile(['sampleneurons_info' mat2str(figurenum) '.pdf']));

%%
% plot the bar plot of recency index and surprise index for one neuron
sampleneuron = Neuronlist_good(Sampleneuron_set(id));

