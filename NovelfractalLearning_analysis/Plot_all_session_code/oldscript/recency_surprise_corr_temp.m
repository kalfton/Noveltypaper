% a tempory script to calculate the correlations of sensory surprise and
% recency in both monkeys
Monkey = 'Lemmy'; %Lemmy, combine

if strcmpi(Monkey, 'Slayer')
    includedneuron = strcmpi({Neuronlist_good(:).monkey}, Monkey);
elseif strcmpi(Monkey, 'Lemmy')
    includedneuron = strcmpi({Neuronlist_good(:).monkey}, Monkey);
elseif strcmpi(Monkey, 'Combine')
    includedneuron = true(size(Neuronlist_good));
else
    error('Wrong Monkey input');
end

clear indices
%lets grab the key indices for hypotheses testing
indices.pred_nov_vs_fam = [Neuronlist_good(includedneuron).pred_nov_vs_fam]';
indices.pred_vs_unpred_fam=[Neuronlist_good(includedneuron).pred_vs_unpred_fam]';
indices.violation_ind=[Neuronlist_good(includedneuron).violation_ind]';
indices.recency_ind=[Neuronlist_good(includedneuron).recency_ind_match_pos]';
indices.uncertaintyindex = [Neuronlist_good(includedneuron).uncertaintyindex_precue_split]';
indices.rewardvalueindex = [Neuronlist_good(includedneuron).rewardvalueindex_precue]';
indices.rewardcuePE = [Neuronlist_good(includedneuron).signedRPE_cue]';
indices.RewInfoAnticipIndex = [Neuronlist_good(includedneuron).RewInfoAnticipIndex_split]';

%lets grab the key indices p values

indices.Ppred_nov_vs_fam = [Neuronlist_good(includedneuron).P_pred_nov_vs_fam]';
indices.Ppred_vs_unpred_fam=[Neuronlist_good(includedneuron).P_pred_vs_unpred_fam_perm]';
indices.Pviolation_ind=[Neuronlist_good(includedneuron).P_violation_ind_perm]';
indices.Precency_ind=[Neuronlist_good(includedneuron).P_recency_ind_match_pos]';
indices.Puncertaintyindex = [Neuronlist_good(includedneuron).uncertaintyindexP_precue_split]';
indices.Prewardvalueindex = [Neuronlist_good(includedneuron).rewardvalueindexP_precue]';
indices.PrewardcuePE = [Neuronlist_good(includedneuron).signedRPE_cue_P]';
indices.PRewInfoAnticipIndex = [Neuronlist_good(includedneuron).RewInfoAnticipIndexP_split]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

heatmapplot_func(indices, plotpath, {['Indices_heatmap_all_session_recency_surprise_' Monkey '.pdf']});
