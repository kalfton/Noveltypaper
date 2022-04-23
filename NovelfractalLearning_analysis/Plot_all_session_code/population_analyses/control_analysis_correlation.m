%% Control analysis for the correlation main bar plots
plotpath = '.\plots';

novelty_indices = {'pred_nov_vs_fam', 'pred_nov_vs_fam_control1', 'pred_nov_vs_fam_control2'};
novelty_indices_P = {'P_pred_nov_vs_fam', 'P_pred_nov_vs_fam_control1', 'P_pred_nov_vs_fam_control2'};
surprise_indices = {'pred_vs_unpred_fam', 'pred_vs_unpred_fam_old', 'pred_vs_unpred_fam_control'};
surprise_indices_P = {'P_pred_vs_unpred_fam_perm', 'P_pred_vs_unpred_fam_perm_old', 'P_pred_vs_unpred_fam_perm_control'};

    
for i_nov = 1:numel(novelty_indices)
    % bar plot analysis
    indices.pred_nov_vs_fam = [Neuronlist_good(:).(novelty_indices{i_nov})]';
    indices.pred_vs_unpred_fam=[Neuronlist_good(:).(surprise_indices{1})]';
    
    indices.Ppred_nov_vs_fam = [Neuronlist_good(:).(novelty_indices_P{i_nov})]';
    indices.Ppred_vs_unpred_fam=[Neuronlist_good(:).(surprise_indices_P{1})]';
    
    plotname = ['Indices_barplot_all_session_' novelty_indices{i_nov} '_' surprise_indices{1} '.pdf'];
    barplot_func(indices, plotpath, plotname);
%     % correlation(error bar) analysis.
%     plotname = {['Indices_heatmap_all_session_' novelty_indices{i_nov} '_' surprise_indices{1} '.pdf']
%         ['Indices_binned_errorbar_all_session' novelty_indices{i_nov} '_' surprise_indices{1} '.pdf']};
%     heatmapplot_func(indices, plotpath, plotname);
    
    

end

% surprise indices
for i_sup = 2:numel(surprise_indices)
    % bar plot analysis
    indices.pred_nov_vs_fam = [Neuronlist_good(:).(novelty_indices{1})]';
    indices.pred_vs_unpred_fam=[Neuronlist_good(:).(surprise_indices{i_sup})]';
    
    indices.Ppred_nov_vs_fam = [Neuronlist_good(:).(novelty_indices_P{1})]';
    indices.Ppred_vs_unpred_fam=[Neuronlist_good(:).(surprise_indices_P{i_sup})]';
    
    plotname = ['Indices_barplot_all_session_' novelty_indices{1} '_' surprise_indices{i_sup} '.pdf'];
    barplot_func(indices, plotpath, plotname);
%     % correlation(error bar) analysis.
%     plotname = {['Indices_heatmap_all_session_' novelty_indices{1} '_' surprise_indices{i_sup} '.pdf']
%         ['Indices_binned_errorbar_all_session' novelty_indices{1} '_' surprise_indices{i_sup} '.pdf']};
%     heatmapplot_func(indices, plotpath, plotname);
    

end

% set the indices back to default:

indices.pred_nov_vs_fam = [Neuronlist_good(:).pred_nov_vs_fam]';
indices.pred_vs_unpred_fam=[Neuronlist_good(:).pred_vs_unpred_fam]';

indices.Ppred_nov_vs_fam = [Neuronlist_good(:).P_pred_nov_vs_fam]';
indices.Ppred_vs_unpred_fam=[Neuronlist_good(:).P_pred_vs_unpred_fam_perm]';


