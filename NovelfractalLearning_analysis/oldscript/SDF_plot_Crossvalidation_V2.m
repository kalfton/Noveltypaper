% SDF plot
figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation

yaxislim = [-0.5,1];
p_thr = StatisticalThreshold;

title_for_plot = {'Novelty', 'Recency', 'pred vs unpred familiar', 'Violation'};
indices_in_struct = {'pred_nov_vs_fam', 'recency_ind_match_pos', 'pred_vs_unpred_fam', 'violation_ind'};
Pvalues_in_struct = {'P_pred_nov_vs_fam', 'P_recency_ind_match_pos', 'P_pred_vs_unpred_fam_perm', 'P_violation_ind_perm'};
SDF_for_plot = {'pred_nov_vs_fam_sdfs','recency_ind_sdfs', 'pred_vs_unpred_fam_sdfs', 'violation_ind_comp_sdf'};
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
        notnan_logic = cellfun(@(x) ~any(isnan(x(:))), {Neuronlist_even(:).(SDF_for_plot{xy})});
        sig_logi = sig_logi & notnan_logic;
        SDF_even_sig = vertcat(Neuronlist_even(sig_logi).(SDF_for_plot{xy}));
        
        P_even = [Neuronlist_even(:).(Pvalues_in_struct{xy})];
        Ind_even = [Neuronlist_even(:).(indices_in_struct{xy})];
        if strcmpi(neurontype{xx}, 'pos')
            sig_logi = P_even<p_thr & Ind_even>0;
        elseif strcmpi(neurontype{xx}, 'neg')
            sig_logi = P_even<p_thr & Ind_even<0;
        end
        notnan_logic = cellfun(@(x) ~any(isnan(x(:))), {Neuronlist_odd(:).(SDF_for_plot{xy})});
        sig_logi = sig_logi & notnan_logic;
        SDF_odd_sig = vertcat(Neuronlist_odd(sig_logi).(SDF_for_plot{xy}));
        
        
        Neurondata = [SDF_odd_sig;SDF_even_sig];
        % substract baseline
        baseline_FR = mean(Neurondata(:,1:100),2);
        baseline_matrix = baseline_FR*ones(1, size(Neurondata,2));
        Neurondata = Neurondata-baseline_matrix;
        
        title([title_for_plot{xy} ' ' neurontype{xx} ' neurons']);
        if size(Neurondata,1)>=4
            x = 1:size(Neurondata,2);
            v = sqrt(size(Neurondata,1)/2);
            plt = shadedErrorBar(x, Neurondata(1:2:end,:), {@nanmean, @(x) nanstd(x)./v}, {'-b', 'LineWidth', 1}, 0);
            set(plt.patch,'FaceAlpha',0.5);
            plt = shadedErrorBar(x, Neurondata(2:2:end,:), {@nanmean, @(x) nanstd(x)./v}, {'-r', 'LineWidth', 1}, 0);
            set(plt.patch,'FaceAlpha',0.5);
        end
        xlim([0,600]);
        ylim(yaxislim);
        y_lim = get(gca,'ylim');
        plot([100,100],y_lim,'k');
        text(200, (0.9*y_lim(1)+0.1*y_lim(2)), sprintf('n = %d', size(Neurondata,1)/2))
        % statistical test
        Mean_FRs = mean(Neurondata(:,101:601),2);
        p = ranksum(Mean_FRs(1:2:end), Mean_FRs(2:2:end));
        text(200, (0.1*y_lim(1)+0.9*y_lim(2)), sprintf('p = %.3d', p));
        
    end
end

%% Preselect the neurons by the responding to novelty

yaxislim = [-0.3,0.4];

for xx = 1:numel(neurontype)
    
    for xy = 1:length(title_for_plot)
        if strcmpi(neurontype{xx}, 'pos')
            nsubplot(169,169, 60+plotplacesety{xy}, plotplacesetx{xy});
        elseif strcmpi(neurontype{xx}, 'neg')
            nsubplot(169,169, 90+plotplacesety{xy}, plotplacesetx{xy});
        else
            error('Invalid neurontype value');
        end
        
        P_odd = [Neuronlist_odd(:).(Pvalues_in_struct{1})];
        Ind_odd = [Neuronlist_odd(:).(indices_in_struct{1})];
        if strcmpi(neurontype{xx}, 'pos')
            sig_logi = P_odd<p_thr & Ind_odd>0;
        elseif strcmpi(neurontype{xx}, 'neg')
            sig_logi = P_odd<p_thr & Ind_odd<0;
        end
        notnan_logic = cellfun(@(x) ~any(isnan(x(:))), {Neuronlist_even(:).(SDF_for_plot{xy})});
        sig_logi = sig_logi & notnan_logic;
        SDF_even_sig = vertcat(Neuronlist_even(sig_logi).(SDF_for_plot{xy}));
        
        P_even = [Neuronlist_even(:).(Pvalues_in_struct{1})];
        Ind_even = [Neuronlist_even(:).(indices_in_struct{1})];
        if strcmpi(neurontype{xx}, 'pos')
            sig_logi = P_even<p_thr & Ind_even>0;
        elseif strcmpi(neurontype{xx}, 'neg')
            sig_logi = P_even<p_thr & Ind_even<0;
        end
        % excude nan sdfs in inhibited neuron
        notnan_logic = cellfun(@(x) ~any(isnan(x(:))), {Neuronlist_odd(:).(SDF_for_plot{xy})});
        sig_logi = sig_logi & notnan_logic;
        SDF_odd_sig = vertcat(Neuronlist_odd(sig_logi).(SDF_for_plot{xy}));

        Neurondata = [SDF_odd_sig;SDF_even_sig];
        % substract baseline
        baseline_FR = mean(Neurondata(:,1:100),2);
        baseline_matrix = baseline_FR*ones(1, size(Neurondata,2));
        Neurondata = Neurondata-baseline_matrix;
        
        title([title_for_plot{xy} ' ' neurontype{xx} ' neurons']);
        if size(Neurondata,1)>=4
            x = 1:size(Neurondata,2);
            v = sqrt(size(Neurondata,1)/2);
            plt = shadedErrorBar(x, Neurondata(1:2:end,:), {@nanmean, @(x) nanstd(x)./v}, {'-b', 'LineWidth', 1}, 0);
            set(plt.patch,'FaceAlpha',0.5);
            plt = shadedErrorBar(x, Neurondata(2:2:end,:), {@nanmean, @(x) nanstd(x)./v}, {'-r', 'LineWidth', 1}, 0);
            set(plt.patch,'FaceAlpha',0.5);
        end
        xlim([0,600]);
        ylim(yaxislim);
        y_lim = get(gca,'ylim');
        plot([100,100],y_lim,'k');
        text(200, (0.9*y_lim(1)+0.1*y_lim(2)), sprintf('n = %d', size(Neurondata,1)/2))
        
        % statistical test
        Mean_FRs = mean(Neurondata(:,101:601),2);
        p = ranksum(Mean_FRs(1:2:end), Mean_FRs(2:2:end));
        text(200, (0.1*y_lim(1)+0.9*y_lim(2)), sprintf('p = %.3d', p));
    end
end

print(gcf,'-dpdf', '-painters',fullfile(plotpath,['SDF_plot_V2_new_crossvalidate.pdf']));