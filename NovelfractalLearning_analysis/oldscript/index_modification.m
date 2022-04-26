function tempneuronlist = index_modification(tempneuronlist, infotempneuronlist)

% add some reward indices into the neuron struct, and rescale some indices to
% [-1, 1].

for xxx = 1:length(tempneuronlist)
    if exist('infotempneuronlist','var')
        % find the spike name in infonew is the same as that in NFL
        matchInd = find(contains({infotempneuronlist(:).name},tempneuronlist(xxx).name));
    else
        matchInd = [];
    end
    
    % If the for the same neuron mua name, the depth infomation is
    % different between the two task( Neuronlist, exclude the neuron)
    if ~isempty(matchInd)
        if infotempneuronlist(matchInd).electrodeDepth ~= tempneuronlist(xxx).electrodeDepth || infotempneuronlist(matchInd).electrodeID ~= tempneuronlist(xxx).electrodeID
            infotempneuronlist(matchInd).exclude = 1;
            infotempneuronlist(matchInd).excludemua = 0;
            tempneuronlist(xxx).exclude = 1;
            tempneuronlist(xxx).excludemua = 0;
        end
    end
    
    if exist('infotempneuronlist','var') && isempty(matchInd)
        error('Neurons should match each other in the two files');
    end
    
    if ~isempty(matchInd) %% && contains(tempneuronlist(xxx).name,'MUA')
        tempneuronlist(xxx).uncertaintyindex = infotempneuronlist(matchInd).uncertaintyindex;
        tempneuronlist(xxx).uncertaintyindexP = infotempneuronlist(matchInd).uncertaintyindexP;
        tempneuronlist(xxx).uncertaintyindex_precue = infotempneuronlist(matchInd).uncertaintyindex_precue;
        tempneuronlist(xxx).uncertaintyindexP_precue = infotempneuronlist(matchInd).uncertaintyindexP_precue;
        tempneuronlist(xxx).uncertaintyindex_preoutcome = infotempneuronlist(matchInd).uncertaintyindex_preoutcome;
        tempneuronlist(xxx).uncertaintyindexP_preoutcome = infotempneuronlist(matchInd).uncertaintyindexP_preoutcome;
        tempneuronlist(xxx).rewardvalueindex = infotempneuronlist(matchInd).rewardvalueindex;
        tempneuronlist(xxx).rewardvalueindexP = infotempneuronlist(matchInd).rewardvalueindexP;
        tempneuronlist(xxx).rewardoutcomePE = infotempneuronlist(matchInd).rewardoutcomePE;
        tempneuronlist(xxx).rewardoutcomePEP = infotempneuronlist(matchInd).rewardoutcomePEP;
        tempneuronlist(xxx).rewardcuePE = infotempneuronlist(matchInd).rewardcuePE;
        tempneuronlist(xxx).rewardcuePEP = infotempneuronlist(matchInd).rewardcuePEP;
        tempneuronlist(xxx).outcometrigSurprise100v50del = infotempneuronlist(matchInd).outcometrigSurprise100v50del;
        tempneuronlist(xxx).outcometrigSurprise100v50delP = infotempneuronlist(matchInd).outcometrigSurprise100v50delP;
        tempneuronlist(xxx).outcometrigSurprise0v50ndel = infotempneuronlist(matchInd).outcometrigSurprise0v50ndel;
        tempneuronlist(xxx).outcometrigSurprise0v50ndelP = infotempneuronlist(matchInd).outcometrigSurprise0v50ndelP;
        
        tempneuronlist(xxx).rewardvalueindex_precue = infotempneuronlist(xxx).rewardvalueindex_precue;
        tempneuronlist(xxx).rewardvalueindexP_precue = infotempneuronlist(xxx).rewardvalueindexP_precue;
        tempneuronlist(xxx).uncertaintyindex_precue_split = infotempneuronlist(xxx).uncertaintyindex_precue_split;
        tempneuronlist(xxx).uncertaintyindexP_precue_split = infotempneuronlist(xxx).uncertaintyindexP_precue_split;
        tempneuronlist(xxx).signedRPE_cuerm_precueroc = infotempneuronlist(xxx).signedRPE_cuerm_precueroc;
        tempneuronlist(xxx).signedRPE_cuerm_precueroc_P = infotempneuronlist(xxx).signedRPE_cuerm_precueroc_P;
        tempneuronlist(xxx).signedRPE_cue = infotempneuronlist(xxx).signedRPE_cue;
        tempneuronlist(xxx).signedRPE_cue_P = infotempneuronlist(xxx).signedRPE_cue_P;
        tempneuronlist(xxx).RewInfoAnticipIndex_split = infotempneuronlist(xxx).RewInfoAnticipIndex_split;
        tempneuronlist(xxx).RewInfoAnticipIndexP_split = infotempneuronlist(xxx).RewInfoAnticipIndexP_split;
        
        tempneuronlist(xxx).RewInfoAnticipIndex = infotempneuronlist(matchInd).RewInfoAnticipIndex;
        tempneuronlist(xxx).RewInfoAnticipIndexP = infotempneuronlist(matchInd).RewInfoAnticipIndexP;
        tempneuronlist(xxx).RewUncOutcomeAnticipIndex = infotempneuronlist(matchInd).RewUncOutcomeAnticipIndex;
        tempneuronlist(xxx).RewUncOutcomeAnticipIndexP = infotempneuronlist(matchInd).RewUncOutcomeAnticipIndexP;
        
        %infonew sdfs
        tempneuronlist(xxx).SDF_rew100ni = infotempneuronlist(matchInd).SDF_rew100ni;
        tempneuronlist(xxx).SDF_rew50ni_del = infotempneuronlist(matchInd).SDF_rew50ni_del;
        tempneuronlist(xxx).SDF_rew50ni_ndel = infotempneuronlist(matchInd).SDF_rew50ni_ndel;
        tempneuronlist(xxx).SDF_rew0ni = infotempneuronlist(matchInd).SDF_rew0ni;
        
        tempneuronlist(xxx).SDF_rew100i = infotempneuronlist(matchInd).SDF_rew100i;
        tempneuronlist(xxx).SDF_rew50i_del = infotempneuronlist(matchInd).SDF_rew50i_del;
        tempneuronlist(xxx).SDF_rew50i_ndel = infotempneuronlist(matchInd).SDF_rew50i_ndel;
        tempneuronlist(xxx).SDF_rew0i = infotempneuronlist(matchInd).SDF_rew0i;
        
        tempneuronlist(xxx).zscorestd_info = infotempneuronlist(matchInd).zscorestd;
        tempneuronlist(xxx).zscoremean_info = infotempneuronlist(matchInd).zscoremean;
    else
        tempneuronlist(xxx).uncertaintyindex = NaN;
        tempneuronlist(xxx).uncertaintyindexP = NaN;
        tempneuronlist(xxx).uncertaintyindex_precue = NaN;
        tempneuronlist(xxx).uncertaintyindexP_precue = NaN;
        tempneuronlist(xxx).uncertaintyindex_preoutcome =NaN;
        tempneuronlist(xxx).uncertaintyindexP_preoutcome = NaN;
        tempneuronlist(xxx).rewardvalueindex = NaN;
        tempneuronlist(xxx).rewardvalueindexP = NaN;
        tempneuronlist(xxx).rewardoutcomePE = NaN;
        tempneuronlist(xxx).rewardoutcomePEP = NaN;
        tempneuronlist(xxx).rewardcuePE = NaN;
        tempneuronlist(xxx).rewardcuePEP = NaN;
        tempneuronlist(xxx).outcometrigSurprise100v50del = NaN;
        tempneuronlist(xxx).outcometrigSurprise100v50delP = NaN;
        tempneuronlist(xxx).outcometrigSurprise0v50ndel =NaN;
        tempneuronlist(xxx).outcometrigSurprise0v50ndelP = NaN;
        
        tempneuronlist(xxx).rewardvalueindex_precue = NaN;
        tempneuronlist(xxx).rewardvalueindexP_precue = NaN;
        tempneuronlist(xxx).uncertaintyindex_precue_split = NaN;
        tempneuronlist(xxx).uncertaintyindexP_precue_split = NaN;
        tempneuronlist(xxx).signedRPE_cuerm_precueroc = NaN;
        tempneuronlist(xxx).signedRPE_cuerm_precueroc_P = NaN;
        tempneuronlist(xxx).signedRPE_cue = NaN;
        tempneuronlist(xxx).signedRPE_cue_P = NaN;
        tempneuronlist(xxx).RewInfoAnticipIndex_split = NaN;
        tempneuronlist(xxx).RewInfoAnticipIndexP_split = NaN;
        
        tempneuronlist(xxx).RewInfoAnticipIndex = NaN;
        tempneuronlist(xxx).RewInfoAnticipIndexP = NaN;
        tempneuronlist(xxx).RewUncOutcomeAnticipIndex = NaN;
        tempneuronlist(xxx).RewUncOutcomeAnticipIndexP = NaN;
        
        %infonew sdfs
        tempneuronlist(xxx).SDF_rew100ni = NaN;
        tempneuronlist(xxx).SDF_rew50ni_del = NaN;
        tempneuronlist(xxx).SDF_rew50ni_ndel = NaN;
        tempneuronlist(xxx).SDF_rew0ni = NaN;
        
        tempneuronlist(xxx).SDF_rew100i = NaN;
        tempneuronlist(xxx).SDF_rew50i_del = NaN;
        tempneuronlist(xxx).SDF_rew50i_ndel = NaN;
        tempneuronlist(xxx).SDF_rew0i = NaN;
        
        tempneuronlist(xxx).zscorestd_info = NaN;
        tempneuronlist(xxx).zscoremean_info = NaN;
    end
    
    %% rescale all necessary indices
    tempneuronlist(xxx).pred_nov_vs_fam = tempneuronlist(xxx).pred_nov_vs_fam*2-1;
    tempneuronlist(xxx).pred_vs_unpred_fam = tempneuronlist(xxx).pred_vs_unpred_fam*2-1;
    tempneuronlist(xxx).violation_ind = tempneuronlist(xxx).violation_ind*2-1;
    tempneuronlist(xxx).recency_ind_match_pos = tempneuronlist(xxx).recency_ind_match_pos*2-1;
    tempneuronlist(xxx).uncertaintyindex_precue_split = tempneuronlist(xxx).uncertaintyindex_precue_split*2;
    tempneuronlist(xxx).rewardvalueindex_precue = tempneuronlist(xxx).rewardvalueindex_precue*2;
    tempneuronlist(xxx).signedRPE_cue = tempneuronlist(xxx).signedRPE_cue*2;
    tempneuronlist(xxx).RewInfoAnticipIndex_split = tempneuronlist(xxx).RewInfoAnticipIndex_split;
    
    tempneuronlist(xxx).pred_nov_vs_fam_control1 = tempneuronlist(xxx).pred_nov_vs_fam_control1*2-1;
    tempneuronlist(xxx).pred_nov_vs_fam_control2 = tempneuronlist(xxx).pred_nov_vs_fam_control2*2-1;
    tempneuronlist(xxx).pred_vs_unpred_fam_old = tempneuronlist(xxx).pred_vs_unpred_fam_old*2-1;
    tempneuronlist(xxx).pred_vs_unpred_fam_control = tempneuronlist(xxx).pred_vs_unpred_fam_control*2-1;
    
    
end

end