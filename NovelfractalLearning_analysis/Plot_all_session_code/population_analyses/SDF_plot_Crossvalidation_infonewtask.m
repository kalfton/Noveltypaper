% SDF plot
figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation

yaxislim = [-0.5,1];
p_thr = StatisticalThreshold;

title_for_plot = {'Reward value', 'Info anticipation info', 'Info anticipation noinfo'};
indices_in_struct = {'rewardvalueindex_precue', 'RewInfoAnticipIndex_split', 'RewInfoAnticipIndex_split'};
Pvalues_in_struct = {'rewardvalueindexP_precue', 'RewInfoAnticipIndexP_split', 'RewInfoAnticipIndexP_split'};


sdfname = {'SDF_rew100ni'
'SDF_rew50ni_del'
'SDF_rew50ni_ndel'
'SDF_rew0ni'
'SDF_rew100i'
'SDF_rew50i_del'
'SDF_rew50i_ndel'
'SDF_rew0i'};
sdf_ind_pair = { {[8,4],[5,1]}, {[8,5],[6,7]}, {[1,4],[2,3]} };

sdf_legend = {{'noreward info', 'reward info'},{'certain info', 'uncertain info'}, {'certain noinfo', 'uncertain noinfo'}};

plotplacesetx = {40:60,70:90,100:120,130:150};
plotplacesety = {5:25, 5:25, 5:25, 5:25 };

neurontype = {'pos', 'neg'};
%% 
for xx = 1:numel(neurontype)
    
    for xy = 1:length(title_for_plot)
        if strcmpi(neurontype{xx}, 'pos')
            nsubplot(169,169, plotplacesety{xy}, plotplacesetx{xy});
        elseif strcmpi(neurontype{xx}, 'neg')
            nsubplot(169,169, 30+plotplacesety{xy}, plotplacesetx{xy});
        else
            error('Invalid neurontype value');
        end
        
        P_odd = [Neuronlist_odd(:).(Pvalues_in_struct{xy})];
        Ind_odd = [Neuronlist_odd(:).(indices_in_struct{xy})];
        if strcmpi(neurontype{xx}, 'pos')
            sig_logi = P_odd<p_thr & Ind_odd>0;
        elseif strcmpi(neurontype{xx}, 'neg')
            sig_logi = P_odd<p_thr & Ind_odd<0;
        end
        %Extract the sdf of selective neurons here.
        sdfs_1 = [];
        sdfs_2 = [];
        
        % Kick out the sdfs that are bad, e.g. nan's or just one number
        notnan_logic = true(size(Neuronlist_even));
        for iii = 1:numel(sdfname)
            notnan_logic = notnan_logic & cellfun(@(x) ~any(isnan(x(:))) & numel(x)>1, {Neuronlist_even(:).(sdfname{iii})});
        end
        sig_logi = sig_logi & notnan_logic;
        
        current_pair = sdf_ind_pair{xy};
        indset1 = current_pair{1};
        indset2 = current_pair{2};
        
        mean_info = vertcat(Neuronlist_even(sig_logi).zscoremean_info);
        std_info = vertcat(Neuronlist_even(sig_logi).zscorestd_info);
        mean_NFL = vertcat(Neuronlist_even(sig_logi).zscoremean);
        std_NFL = vertcat(Neuronlist_even(sig_logi).zscorestd);
        
        for iii = 1:numel(indset1)
            tempsdfs = vertcat(Neuronlist_even(sig_logi).(sdfname{indset1(iii)}));
            sdfs_new = (tempsdfs.*std_info+mean_info-mean_NFL)./std_NFL;
            sdfs_1 = [sdfs_1; tempsdfs];
        end
        for iii = 1:numel(indset2)
            tempsdfs = vertcat(Neuronlist_even(sig_logi).(sdfname{indset2(iii)}));
            sdfs_new = (tempsdfs.*std_info+mean_info-mean_NFL)./std_NFL;
            sdfs_2 = [sdfs_2; tempsdfs];
        end
        n_neuron = sum(sig_logi);
        
        
        P_even = [Neuronlist_even(:).(Pvalues_in_struct{xy})];
        Ind_even = [Neuronlist_even(:).(indices_in_struct{xy})];
        if strcmpi(neurontype{xx}, 'pos')
            sig_logi = P_even<p_thr & Ind_even>0;
        elseif strcmpi(neurontype{xx}, 'neg')
            sig_logi = P_even<p_thr & Ind_even<0;
        end
        
        % Extract the sdf of selective neurons here.
        
        % Kick out the sdfs that are bad, e.g. nan's or just one number
        notnan_logic = true(size(Neuronlist_odd));
        for iii = 1:numel(sdfname)
            notnan_logic = notnan_logic & cellfun(@(x) ~any(isnan(x(:))) & numel(x)>1, {Neuronlist_odd(:).(sdfname{iii})});
        end
        sig_logi = sig_logi & notnan_logic;
        
        current_pair = sdf_ind_pair{xy};
        indset1 = current_pair{1};
        indset2 = current_pair{2};
        
        mean_info = vertcat(Neuronlist_odd(sig_logi).zscoremean_info);
        std_info = vertcat(Neuronlist_odd(sig_logi).zscorestd_info);
        mean_NFL = vertcat(Neuronlist_odd(sig_logi).zscoremean);
        std_NFL = vertcat(Neuronlist_odd(sig_logi).zscorestd);
        
        for iii = 1:numel(indset1)
            tempsdfs = vertcat(Neuronlist_odd(sig_logi).(sdfname{indset1(iii)}));
            sdfs_new = (tempsdfs.*std_info+mean_info-mean_NFL)./std_NFL;
            sdfs_1 = [sdfs_1; tempsdfs];
        end
        for iii = 1:numel(indset2)
            tempsdfs = vertcat(Neuronlist_odd(sig_logi).(sdfname{indset2(iii)}));
            sdfs_new = (tempsdfs.*std_info+mean_info-mean_NFL)./std_NFL;
            sdfs_2 = [sdfs_2; tempsdfs];
        end
        n_neuron = n_neuron+sum(sig_logi);
        
        
        %Neurondata = [SDF_odd_sig;SDF_even_sig];
        title([title_for_plot{xy} ' ' neurontype{xx} ' neurons']);
        x = 1:size(sdfs_1,2);
        v = sqrt(size(sdfs_1,1)/2); % may need to change it (/2) b.c. it contains 2sdf in each file
        plt1 = shadedErrorBar(x, sdfs_1, {@nanmean, @(x) nanstd(x)./v}, {'-b', 'LineWidth', 1}, 0);
        set(plt1.patch,'FaceAlpha',0.5);
        x = 1:size(sdfs_2,2);
        v = sqrt(size(sdfs_2,1)/2); % may need to change it (/2) b.c. it contains 2sdf in each file
        plt2 = shadedErrorBar(x, sdfs_2, {@nanmean, @(x) nanstd(x)./v}, {'-r', 'LineWidth', 1}, 0);
        set(plt2.patch,'FaceAlpha',0.5);
        
        xlim([0,2000]);
        ylim(yaxislim);
        y_lim = get(gca,'ylim');
        plot([1000,1000],y_lim,'k');
        text(200, (0.9*y_lim(1)+0.1*y_lim(2)), sprintf('n = %d', n_neuron)) %Note this may need to change
        legend([plt1.mainLine, plt2.mainLine], sdf_legend{xy}, 'Fontsize', 5);
    end
end

%% Preselect the neurons by the responding to novelty
for xx = 1:numel(neurontype)
    
    for xy = 1:length(title_for_plot)
        if strcmpi(neurontype{xx}, 'pos')
            nsubplot(169,169, 60+plotplacesety{xy}, plotplacesetx{xy});
        elseif strcmpi(neurontype{xx}, 'neg')
            nsubplot(169,169, 90+plotplacesety{xy}, plotplacesetx{xy});
        else
            error('Invalid neurontype value');
        end
        
        P_odd = [Neuronlist_odd(:).P_pred_nov_vs_fam];
        Ind_odd = [Neuronlist_odd(:).pred_nov_vs_fam];
        if strcmpi(neurontype{xx}, 'pos')
            sig_logi = P_odd<p_thr & Ind_odd>0;
        elseif strcmpi(neurontype{xx}, 'neg')
            sig_logi = P_odd<p_thr & Ind_odd<0;
        end
        %Extract the sdf of selective neurons here.
        sdfs_1 = [];
        sdfs_2 = [];
        
        % Kick out the sdfs that are bad, e.g. nan's or just one number
        notnan_logic = true(size(Neuronlist_even));
        for iii = 1:numel(sdfname)
            notnan_logic = notnan_logic & cellfun(@(x) ~any(isnan(x(:))) & numel(x)>1, {Neuronlist_even(:).(sdfname{iii})});
        end
        sig_logi = sig_logi & notnan_logic;
        
        current_pair = sdf_ind_pair{xy};
        indset1 = current_pair{1};
        indset2 = current_pair{2};
        
        mean_info = vertcat(Neuronlist_even(sig_logi).zscoremean_info);
        std_info = vertcat(Neuronlist_even(sig_logi).zscorestd_info);
        mean_NFL = vertcat(Neuronlist_even(sig_logi).zscoremean);
        std_NFL = vertcat(Neuronlist_even(sig_logi).zscorestd);
        
        for iii = 1:numel(indset1)
            tempsdfs = vertcat(Neuronlist_even(sig_logi).(sdfname{indset1(iii)}));
            sdfs_new = (tempsdfs.*std_info+mean_info-mean_NFL)./std_NFL;
            sdfs_1 = [sdfs_1; tempsdfs];
        end
        for iii = 1:numel(indset2)
            tempsdfs = vertcat(Neuronlist_even(sig_logi).(sdfname{indset2(iii)}));
            sdfs_new = (tempsdfs.*std_info+mean_info-mean_NFL)./std_NFL;
            sdfs_2 = [sdfs_2; tempsdfs];
        end
        n_neuron = sum(sig_logi);
%         %%%
%        
%         notnan_logic = cellfun(@(x) ~any(isnan(x(:))), {Neuronlist_even(:).(SDF_for_plot{xy})});
%         sig_logi = sig_logi & notnan_logic;
%         SDF_even_sig = vertcat(Neuronlist_even(sig_logi).(SDF_for_plot{xy}));
        
        P_even = [Neuronlist_even(:).P_pred_nov_vs_fam];
        Ind_even = [Neuronlist_even(:).pred_nov_vs_fam];
        if strcmpi(neurontype{xx}, 'pos')
            sig_logi = P_even<p_thr & Ind_even>0;
        elseif strcmpi(neurontype{xx}, 'neg')
            sig_logi = P_even<p_thr & Ind_even<0;
        end
        
        % Extract the sdf of selective neurons here.
        
        % Kick out the sdfs that are bad, e.g. nan's or just one number
        notnan_logic = true(size(Neuronlist_odd));
        for iii = 1:numel(sdfname)
            notnan_logic = notnan_logic & cellfun(@(x) ~any(isnan(x(:))) & numel(x)>1, {Neuronlist_odd(:).(sdfname{iii})});
        end
        sig_logi = sig_logi & notnan_logic;
        
        current_pair = sdf_ind_pair{xy};
        indset1 = current_pair{1};
        indset2 = current_pair{2};
        
        mean_info = vertcat(Neuronlist_odd(sig_logi).zscoremean_info);
        std_info = vertcat(Neuronlist_odd(sig_logi).zscorestd_info);
        mean_NFL = vertcat(Neuronlist_odd(sig_logi).zscoremean);
        std_NFL = vertcat(Neuronlist_odd(sig_logi).zscorestd);
        
        for iii = 1:numel(indset1)
            tempsdfs = vertcat(Neuronlist_odd(sig_logi).(sdfname{indset1(iii)}));
            sdfs_new = (tempsdfs.*std_info+mean_info-mean_NFL)./std_NFL;
            sdfs_1 = [sdfs_1; tempsdfs];
        end
        for iii = 1:numel(indset2)
            tempsdfs = vertcat(Neuronlist_odd(sig_logi).(sdfname{indset2(iii)}));
            sdfs_new = (tempsdfs.*std_info+mean_info-mean_NFL)./std_NFL;
            sdfs_2 = [sdfs_2; tempsdfs];
        end
        n_neuron = n_neuron+sum(sig_logi);
        
        
        %Neurondata = [SDF_odd_sig;SDF_even_sig];
        title([title_for_plot{xy} ' ' neurontype{xx} ' neurons']);
        x = 1:size(sdfs_1,2);
        v = sqrt(size(sdfs_1,1)/2); % may need to change it (/2) b.c. it contains 2sdf in each file
        plt1 = shadedErrorBar(x, sdfs_1, {@nanmean, @(x) nanstd(x)./v}, {'-b', 'LineWidth', 1}, 0);
        set(plt1.patch,'FaceAlpha',0.5);
        x = 1:size(sdfs_2,2);
        v = sqrt(size(sdfs_2,1)/2); % may need to change it (/2) b.c. it contains 2sdf in each file
        plt2 = shadedErrorBar(x, sdfs_2, {@nanmean, @(x) nanstd(x)./v}, {'-r', 'LineWidth', 1}, 0);
        set(plt2.patch,'FaceAlpha',0.5);
        
        xlim([0,2000]);
        ylim(yaxislim);
        y_lim = get(gca,'ylim');
        plot([1000,1000],y_lim,'k');
        text(200, (0.9*y_lim(1)+0.1*y_lim(2)), sprintf('n = %d', n_neuron)) %Note this may need to change
        legend([plt1.mainLine, plt2.mainLine], sdf_legend{xy}, 'Fontsize', 5);
    end
end

print(gcf,'-dpdf', '-painters',fullfile(plotpath,['SDF_plot_info_new_crossvalidate.pdf']));

