%%
%% Second plot to check the q value/ distribution of p value for info task 
figure;
set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation


Neuronlist_infonew_good = Neuronlist_infonew_all([Neuronlist_infonew_all(:).exclude]==0 | [Neuronlist_infonew_all(:).excludemua]==1);
%% plot P value;

Pvaluenames = {'rewardvalueindexP_precue','RewInfoAnticipIndexP_split'};
Pvaluelabel = {'Reward value','Infoanticip'};
Qvaluenames = {'Q_rewardvalueindex_info','Q_RewInfoAnticipIndex_info'};

    
plotplacesetx = {101:120, 126:145, 151:170, 176:195};
plotplacesety = {1:30,1:30,1:30, 1:30};
plotspace = 50;

for vw = 1:length(Pvaluenames)

Pvalue = [Neuronlist_infonew_good(:).(Pvaluenames{vw})]; 
%Too far indices: P_rocObjSparseness

%Good indices: P_pred_nov_vs_fam, P_position_eff_anova, P_unp_nov_vs_fam,
%P_violation_ind_perm, P_reward_short_vs_long

%Between: P_pred_vs_unpred_fam_perm

%Bad indices: P_recency_ind, P_object_select_anova

%%% if all p values is nan, skip the Q value calculation, which is the case
%%% when exlude mua and try to use reward related p value
if sum(~isnan(Pvalue))<=1
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


print(gcf,'-dpdf', '-painters',fullfile(plotpath,['Qvalue_info_all_session' '.pdf']));




