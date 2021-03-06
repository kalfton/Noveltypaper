%shuffling_num = 10000;
run_longtimebootstrapping = 0;
% 
%Monkey S
fractalIDset_S = {6300:6303, 7999, [7410:7411,7420:7421], [7412,7422,7413,7423, 7414,7424, 7415,7425]};
%Monkey L
fractalIDset_L = {6300:6307, 7999, [7300:7307], [7410:7411,7420:7421]};

% choose the neurons in the sessions which has multiday fractal
logical_multiday = cellfun(@(x) ~isempty(x.('FR7410')) | ~isempty(x.('FR7411')), {Neuronlist_all(:).learning})';
logical_multiday = logical_multiday & cellfun(@(x) (numel(x.('learningdate'))==5 || numel(x.('learningdate'))==1 && x.('learningdate')>1), {Neuronlist_all(:).learning})';
logical_multiday = logical_multiday;

%%
%% repeating novel fractals plot by presentation number.

logical_for_neuronlist = {[Neuronlist_all(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_all(:).pred_nov_vs_fam]'>0 & logical_multiday
     };
    %[Neuronlist_all(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_all(:).pred_nov_vs_fam]'<0  & logical_multiday


Select_criteria = {'Nov excited', 'Nov Inhibited'};

Regression_x_names = {'Regression_FL_x', 'Regression_NL_x', 'Regression_L1_x', 'Regression_L2_x', 'Regression_L3_x', 'Regression_L4_x', 'Regression_L5_x'};
All_FR_names = {'All_FR_FL', 'All_FR_NL', 'All_FR_L1', 'All_FR_L2', 'All_FR_L3', 'All_FR_L4', 'All_FR_L5'};
exp_paras_names = {'exp_paras_FL', 'exp_paras_NL', 'exp_paras_L1', 'exp_paras_L2', 'exp_paras_L3', 'exp_paras_L4', 'exp_paras_L5'};

for xy = 1:numel(logical_for_neuronlist) % novelty excited, inhibited
    

Neuronlist_learning = Neuronlist_all(logical_for_neuronlist{xy});
num = sum(logical_for_neuronlist{xy});
figure;
set(gcf,'Position',[100 100 1100 800],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
nsubplot(3,3,1,1);
hold on;
% initialize
fractalIDset = fractalIDset_S; % default
for xxx = 1:length(fractalIDset)
    eval([Regression_x_names{xxx} ' = [];']);
    eval([All_FR_names{xxx} ' = [];']);
end


for xxx = 1:length(fractalIDset)
    Regression_x = [];
    All_FR = [];
    for iii = 1:length(Neuronlist_learning)
        if strcmpi(Neuronlist_learning(iii).monkeyName, 'L')
            fractalIDset = fractalIDset_L;
        elseif strcmpi(Neuronlist_learning(iii).monkeyName, 'S')
            fractalIDset = fractalIDset_S;
        end
        
        for ii = 1:length(fractalIDset{xxx})
            if isfield(Neuronlist_learning(iii).learning, ['FR' mat2str(fractalIDset{xxx}(ii))])
                Regression_x = [Regression_x; (1:length(Neuronlist_learning(iii).learning.(['FR' mat2str(fractalIDset{xxx}(ii))])))'];
                All_FR = [All_FR; Neuronlist_learning(iii).learning.(['FR' mat2str(fractalIDset{xxx}(ii))])];
            end
        end
    end
    eval([Regression_x_names{xxx} '= Regression_x;']);
    eval([All_FR_names{xxx} '= All_FR;']);
end


nsubplot(3,3,1,1);
meanFR_FL = zeros(1,max(Regression_L1_x));
meanFR_NL = zeros(1,max(Regression_L1_x));
meanFR_L1 = zeros(1,max(Regression_L1_x));
meanFR_L2 = zeros(1,max(Regression_L1_x));
try
    meanFR_L3 = zeros(1,max(Regression_L1_x));
    meanFR_L4 = zeros(1,max(Regression_L1_x));
    meanFR_L5 = zeros(1,max(Regression_L1_x));
catch
end

for ii = 1:max(Regression_L1_x)
    meanFR_FL(ii) = mean(All_FR_FL(Regression_FL_x==ii));
    meanFR_NL(ii) = mean(All_FR_NL(Regression_NL_x==ii));
    meanFR_L1(ii) = mean(All_FR_L1(Regression_L1_x==ii));
    meanFR_L2(ii) = mean(All_FR_L2(Regression_L2_x==ii));
    try
        meanFR_L3(ii) = mean(All_FR_L3(Regression_L3_x==ii));
        meanFR_L4(ii) = mean(All_FR_L4(Regression_L4_x==ii));
        meanFR_L5(ii) = mean(All_FR_L5(Regression_L5_x==ii));
    catch
    end
end
stdFR_FL = zeros(1,max(Regression_L1_x));
stdFR_NL = zeros(1,max(Regression_L1_x));
stdFR_L1 = zeros(1,max(Regression_L1_x));
stdFR_L2 = zeros(1,max(Regression_L1_x));
for ii = 1:max(Regression_L1_x)
    stdFR_FL(ii) = std(All_FR_FL(Regression_FL_x==ii))/sqrt(sum(Regression_FL_x==ii));
    stdFR_NL(ii) = std(All_FR_NL(Regression_NL_x==ii))/sqrt(sum(Regression_FL_x==ii));
    stdFR_L1(ii) = std(All_FR_L1(Regression_L1_x==ii))/sqrt(sum(Regression_FL_x==ii));
    stdFR_L2(ii) = std(All_FR_L2(Regression_L2_x==ii))/sqrt(sum(Regression_FL_x==ii));
end

plot_x = 1:max(Regression_L1_x);
hold on;
plot(nan, nan,'r');
plot(nan, nan,'k');
plot(nan, nan,'b');
plot(nan, nan,'color',[0,0.3,0.7]);
shadedErrorBar(plot_x,meanFR_NL,stdFR_NL,{'-r', 'LineWidth', 0.5}, 0)
shadedErrorBar(plot_x,meanFR_FL,stdFR_FL,{'-k', 'LineWidth', 0.5}, 0)
shadedErrorBar(plot_x,meanFR_L1,stdFR_L1,{'-b', 'LineWidth', 0.5}, 0)
shadedErrorBar(plot_x,meanFR_L2,stdFR_L2,{'color',[0,0.3,0.7], 'LineWidth', 0.5}, 0)

xlim([0,20]);
ylim([-0.5,1]);
lgd = legend('novel', 'familiar', 'learning 1 day', 'learning n>=2 day');
lgd.Position(1) = lgd.Position(1)+0.2;
xlabel('Appearance times')
ylabel('Z-scored firing rate');
title([Select_criteria{xy} ', n=' mat2str(num)]) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Normalizing firing rate and measuring learning speed
normalized_meanFR_L1 = (meanFR_L1-meanFR_FL)./(meanFR_NL-meanFR_FL);
normalized_meanFR_L2 = (meanFR_L2-meanFR_FL)./(meanFR_NL-meanFR_FL);

learning_rate_L1 = -(normalized_meanFR_L1(2:end)-normalized_meanFR_L1(1:end-1))./normalized_meanFR_L1(1:end-1);
learning_rate_L2 = -(normalized_meanFR_L2(2:end)-normalized_meanFR_L2(1:end-1))./normalized_meanFR_L2(1:end-1);

smoothfun = @(x) smooth(x,3,'moving');
learning_rate_L1_s = smoothfun(learning_rate_L1)';%smooth(learning_rate_L1,3,'moving');
learning_rate_L2_s = smoothfun(learning_rate_L2)';%smooth(learning_rate_L2,3,'moving');

if xy>1 || ~run_longtimebootstrapping
    nsubplot(3,3,2,1);
    hold on;
    plot_x = 1:max(Regression_L1_x);
    plot(plot_x, normalized_meanFR_L1,'b');
    plot(plot_x, normalized_meanFR_L2,'color',[0,0.3,0.7]);
    ylabel('normalized Firing rate');
    xlim([0,20]);
    ylim([0,1]);
    
    
    nsubplot(3,3,2,2);
    plot_x = 1:max(Regression_L1_x)-1;
    
    hold on;
    plot(plot_x, learning_rate_L1_s,'b');
    plot(plot_x, learning_rate_L2_s,'color',[0,0.3,0.7]);
    
    ylabel('Learning rate (alpha)')
    xlim([0,20]);
    
%     %gaussian process regression
%     sigma0 = 0.1;
%     kparams0 = [0.05, 0.2];
%     gprMd2 = fitrgp(plot_x(1:20)', learning_rate_L1(1:20),'KernelFunction','squaredexponential',...
%          'KernelParameters', kparams0, 'Sigma',sigma0);
%     
%     ypred2 = predict(gprMd2,plot_x');
%     
%     plot(plot_x, ypred2,'r');
    
else
    %% bootstrapping to get confidence interval: 
    
    pars.shuffling_num = shuffling_num;
    pars.fractalIDset_S = fractalIDset_S;
    pars.fractalIDset_L = fractalIDset_L;
    pars.maxappearance = max(Regression_L1_x);
    pars.smoothfun = smoothfun;
    
    datastruct = Learning_acrossday_bootstrapping(Neuronlist_learning,pars);
    
    normalized_meanFR_L1_bootstrapping = datastruct.normalized_meanFR_L1_bootstrapping;
    normalized_meanFR_L2_bootstrapping = datastruct.normalized_meanFR_L2_bootstrapping;
    learning_rate_L1_s_bootstrapping = datastruct.learning_rate_L1_s_bootstrapping;
    learning_rate_L2_s_bootstrapping = datastruct.learning_rate_L2_s_bootstrapping;
        
    maxappearance = max(Regression_L1_x);

    normalized_meanFR_L1_bootstrapping = sort(normalized_meanFR_L1_bootstrapping,1);
    normalized_meanFR_L2_bootstrapping = sort(normalized_meanFR_L2_bootstrapping,1);
    learning_rate_L1_s_bootstrapping = sort(learning_rate_L1_s_bootstrapping,1);
    learning_rate_L2_s_bootstrapping = sort(learning_rate_L2_s_bootstrapping,1);
    
    upperbound_normalized_meanFR_L1 = normalized_meanFR_L1_bootstrapping(ceil(shuffling_num-shuffling_num*0.025), :);
    lowerbound_normalized_meanFR_L1 = normalized_meanFR_L1_bootstrapping(floor(shuffling_num*0.025), :);
    bound_normalized_meanFR_L1 = [upperbound_normalized_meanFR_L1 - normalized_meanFR_L1;...
        normalized_meanFR_L1 - lowerbound_normalized_meanFR_L1];
    
    upperbound_normalized_meanFR_L2 = normalized_meanFR_L2_bootstrapping(ceil(shuffling_num-shuffling_num*0.025), :);
    lowerbound_normalized_meanFR_L2 = normalized_meanFR_L2_bootstrapping(floor(shuffling_num*0.025), :);
    bound_normalized_meanFR_L2 = [upperbound_normalized_meanFR_L2 - normalized_meanFR_L2;...
        normalized_meanFR_L2 - lowerbound_normalized_meanFR_L2];
    
    upperbound_learning_rate_L1_s = learning_rate_L1_s_bootstrapping(ceil(shuffling_num-shuffling_num*0.025), :);
    lowerbound_learning_rate_L1_s = learning_rate_L1_s_bootstrapping(floor(shuffling_num*0.025), :);
    bound_learning_rate_L1_s = [upperbound_learning_rate_L1_s - learning_rate_L1_s;...
        learning_rate_L1_s - lowerbound_learning_rate_L1_s];
    
    upperbound_learning_rate_L2_s = learning_rate_L2_s_bootstrapping(ceil(shuffling_num-shuffling_num*0.025), :);
    lowerbound_learning_rate_L2_s = learning_rate_L2_s_bootstrapping(floor(shuffling_num*0.025), :);
    bound_learning_rate_L2_s = [upperbound_learning_rate_L2_s - learning_rate_L2_s;...
        learning_rate_L2_s - lowerbound_learning_rate_L2_s];
    
    %%% overwrite the errorbar bound with standard deviation.
    bound_normalized_meanFR_L1 = std(normalized_meanFR_L1_bootstrapping);
    bound_normalized_meanFR_L2 = std(normalized_meanFR_L2_bootstrapping);
    bound_learning_rate_L1_s = std(learning_rate_L1_s_bootstrapping);
    bound_learning_rate_L2_s = std(learning_rate_L2_s_bootstrapping);
    
    
    nsubplot(3,3,2,1);
    hold on;
    plot_x = 1:maxappearance;
    shadedErrorBar(plot_x,normalized_meanFR_L1,bound_normalized_meanFR_L1,{'-b', 'LineWidth', 0.5}, 0);
    shadedErrorBar(plot_x,normalized_meanFR_L2,bound_normalized_meanFR_L2,{'color',[0,0.3,0.7], 'LineWidth', 0.5}, 0);
    ylabel('normalized Firing rate');
    xlim([0,20]);
    ylim([0,1]);
    
    
    nsubplot(3,3,2,2);
    plot_x = 1:maxappearance-1;
    hold on;
    shadedErrorBar(plot_x,learning_rate_L1_s,bound_learning_rate_L1_s,{'-b', 'LineWidth', 0.5}, 0);
    shadedErrorBar(plot_x,learning_rate_L2_s,bound_learning_rate_L2_s,{'color',[0,0.3,0.7], 'LineWidth', 0.5}, 0);
    ylabel('Learning rate (alpha)')
    %legend({'Familiar vs Learning 2 days', 'Learning 2 days vs 1 day', '1 day vs Novel'});
    xlim([0,20]);

end
print(gcf,'-dpdf', '-painters',[plotpath '/Learning_n_day_fitting' Select_criteria{xy} '.pdf']);
end


% function sq_err = exp_fun_fitting(Regression_x, exp_paras, firingrate)
% %%% calculate squared error given a set of parameters.
% % Regression_x here should a column vector of x data, firing rate should be
% % a column of y data, exp_paras is the current parameter(length: 3), and the
% % function returns the error of the fitting function with the actual
% % function.
% sq_err = 0;
% for iii = 1:size(Regression_x, 1)
%     appearingtime = Regression_x(iii,1);
%     yvalue = exp_paras(1)+exp_paras(2)*exp(-exp_paras(3)*appearingtime);
%     sq_err = sq_err+(firingrate(iii)-yvalue)^2;
% end
%     
% end