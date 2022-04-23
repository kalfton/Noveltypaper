figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation

%% plot P value;

Pvaluenames = {'P_pred_nov_vs_fam', 'P_recency_ind_match_pos', 'P_violation_ind_perm','P_pred_vs_unpred_fam_perm','uncertaintyindexP_precue_split','rewardvalueindexP_precue','signedRPE_cuerm_precueroc_P', 'RewInfoAnticipIndexP_split'};
Pvaluelabel = {'Novelty', 'Recency', 'Violation', 'Sensory Surprise', 'Uncertainty', 'Reward value', 'RPE','Infoanticip'};
Qvaluenames = {'Q_pred_nov_vs_fam', 'Q_recency_ind', 'Q_violation_ind_perm','Q_pred_vs_unpred_fam_perm','Q_uncertaintyindex','Q_rewardvalueindex','Q_rewardcuePE','Q_RewInfoAnticipIndex'};

plotplacesetx = {1:20,26:45,51:70,76:95, 101:120, 126:145, 151:170, 176:195};
plotplacesety = {1:30,1:30,1:30,1:30, 1:30,1:30,1:30, 1:30};
plotspace = 50;

for vw = 1:length(Pvaluenames)
try
    if strcmp(trialset,'all')
        Pvalue = [Neuronlist_good(:).(Pvaluenames{vw})];
        var_append = '';
    elseif strcmp(trialset,'odd')
        Pvalue = [Neuronlist_odd(:).(Pvaluenames{vw})];
        var_append = '_odd';
    elseif strcmp(trialset,'even')
        Pvalue = [Neuronlist_even(:).(Pvaluenames{vw})];
        var_append = '_even';
    else
        warning('wrong trialset variable value.');
    end
catch
    continue;
end
%Too far indices: P_rocObjSparseness

%Good indices: P_pred_nov_vs_fam, P_position_eff_anova, P_unp_nov_vs_fam,
%P_violation_ind_perm, P_reward_short_vs_long

%Between: P_pred_vs_unpred_fam_perm

%Bad indices: P_recency_ind, P_object_select_anova


%%% if all p values is nan, skip the Q value calculation, which is the case
%%% when exlude mua and try to use reward related p value
if sum(~isnan(Pvalue))<=1
    eval([Qvaluenames{vw} var_append ' = nan(size(Pvalue));']);
    continue;
end

sampleP_threshold = [0.001:0.001:0.01, 0.01:0.01:0.9];
N_significant = zeros(size(sampleP_threshold));
N_population_thrh = zeros(size(sampleP_threshold));
Qvalue_array = zeros(size(sampleP_threshold));

% Calculating Q value
[fdr,q,priori] = mafdr(Pvalue,'Showplot',false,'Method','polynomial');

eval([Qvaluenames{vw} ' = q;']);
[Pvalue,I] = sort(Pvalue);
q = q(I);
Pvalue(isnan(Pvalue))=[];
q(isnan(q))=[];

for iii = 1:length(sampleP_threshold)
    N_significant(iii) = length(find(Pvalue<=sampleP_threshold(iii)));
    % population significant test, the number of significant neuron is
    % approximately poisson distribution
    lambta = length(Pvalue)*sampleP_threshold(iii);
    N_population_thrh(iii) = find(poisscdf([0:length(Pvalue)],lambta)>=0.95, 1,'first');
    
    try
        Qvalue_array(iii) = q(find(Pvalue<sampleP_threshold(iii),1,'last'));
    catch
        Qvalue_array(iii) = q(1);
    end
end
ratio_significant = N_significant/length(Pvalue);
ratio_population_thrh = N_population_thrh/length(Pvalue);


nsubplot(169,195, plotplacesety{vw}, plotplacesetx{vw});
hold on;
plot(sampleP_threshold, ratio_significant,'o');
plot(sampleP_threshold, sampleP_threshold,'r-');
plot(sampleP_threshold, ratio_population_thrh, 'm--');
text(0.05,0.1, ['n=' mat2str(length(Pvalue))]);
xlim([0,0.1]);
ylim([0,0.15]);
xlabel('p value threshold');
ylabel('No. of significant neurons');
title(Pvaluelabel{vw});

% q(p) plot
nsubplot(169,195, plotplacesety{vw}+plotspace, plotplacesetx{vw});
plot(sampleP_threshold,Qvalue_array);
xlim([0,0.1]);
ylim([0,1]);
ylabel('Q value');
xlabel('pvalue');

% binomial test at threshold 0.01
sig =length(find(Pvalue<=0.01));
p = myBinomTest(sig,length(Pvalue),0.01,'Greater');
title(['P binomial(0.01)= ', mat2str(p,3)]);

if vw == 2
    pause(1);
end
% histogram of pvalue

nsubplot(169,195, plotplacesety{vw}+2*plotspace, plotplacesetx{vw});
sampleP_threshold = [0.01:0.01:1];
histogram(Pvalue,[0,sampleP_threshold]);
title('Histogram of P value');
xlim([-0.01,1]);
hold on;
line([0,1], [1,1]*numel(Pvalue)*0.01*priori,'color', 'r', 'LineWidth',1);
% N = histcounts(Pvalue,[0,sampleP_threshold]);
% title(mat2str(sum(N)))

end


print(gcf,'-dpdf', '-painters',fullfile(plotpath,['Qvalue_all_session' var_append '.pdf']));
%% false positive rate:
% Realization of the algorithm of the paper John & Robert, 2003
% 
% 
% figure;
% sampleP_threshold = [0.01:0.01:0.5];
% histogram(Pvalue,[0,sampleP_threshold]);

%number = histcounts(Pvalue,[0,sampleP_threshold]);

% %%% calculate estimated pi0
% pi0_array = zeros(size(sampleP_threshold));
% for iii = 1:length(sampleP_threshold)
%     pi0_array(iii) = sum(Pvalue>sampleP_threshold(iii))/m/(1-sampleP_threshold(iii));
% end
% %%% fitting
% 
% 
% figure;
% plot(sampleP_threshold,  pi0_array);
% ylim([0.85,1]);
%
% figure; plot([fdr;q]');
% legend('fdr','q');


