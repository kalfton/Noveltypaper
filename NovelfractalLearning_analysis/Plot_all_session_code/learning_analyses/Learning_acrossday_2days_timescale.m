close all;
%shuffling_num = 10000;
run_longtimebootstrapping = 0;
% 
%Slayer
fractalIDset_Slayer = {6300:6303, 7999, [7410:7411,7420:7421], [7412,7422,7413,7423, 7414,7424, 7415,7425]};
logical_Slayer = cellfun(@(x) (strcmpi(x ,'Slayer')), {Neuronlist_good(:).monkeyName})';
%Lemmy
fractalIDset_Lemmy = {6300:6307, 7999, [7300:7307], [7410:7411,7420:7421]};
logical_Lemmy = cellfun(@(x) (strcmpi(x ,'Lemmy')), {Neuronlist_good(:).monkeyName})';

% plot the neuron's learning curve separately for novelty excited,
% inhibited, and other neuron.
% And recency
% choose the neurons in the sessions which has multiday fractal
logical_multiday = cellfun(@(x) ~isempty(x.('FR7410')) | ~isempty(x.('FR7411')), {Neuronlist_good(:).learning})';
logical_multiday = logical_multiday & cellfun(@(x) (numel(x.('learningdate'))==5 || numel(x.('learningdate'))==1 && x.('learningdate')>1), {Neuronlist_good(:).learning})';
logical_multiday = logical_multiday;


variablenames = {'Learning_Familiar', 'Learning_Novel', 'Learning_1day', 'Learning_2day'};
fractaldateset = {nan, nan, nan, nan, nan, nan, nan};

%%
%% appearance by appearance plot and exponential fitting

logical_for_neuronlist = {[Neuronlist_good(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]'>0 & logical_multiday ...
     [Neuronlist_good(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]'<0  & logical_multiday ...
     };


Select_criteria = {'Nov excited', 'Nov Inhibited', 'Fast learning fast forgetting', 'Slow learning slow forgetting', 'No learning no forgetting'...
    'Negative learning negative forgetting'};

Regression_x_names = {'Regression_FL_x', 'Regression_NL_x', 'Regression_L1_x', 'Regression_L2_x', 'Regression_L3_x', 'Regression_L4_x', 'Regression_L5_x'};
All_FR_names = {'All_FR_FL', 'All_FR_NL', 'All_FR_L1', 'All_FR_L2', 'All_FR_L3', 'All_FR_L4', 'All_FR_L5'};
exp_paras_names = {'exp_paras_FL', 'exp_paras_NL', 'exp_paras_L1', 'exp_paras_L2', 'exp_paras_L3', 'exp_paras_L4', 'exp_paras_L5'};

for xy = 1:2 % novelty excited, inhibited
    

Neuronlist_learning = Neuronlist_good(logical_for_neuronlist{xy});
num = sum(logical_for_neuronlist{xy});
figure;
nsubplot(75,75,1:20,1:20);
hold on;
% initialize
fractalIDset = fractalIDset_Slayer; % default
for xxx = 1:length(fractalIDset)
    eval([Regression_x_names{xxx} ' = [];']);
    eval([All_FR_names{xxx} ' = [];']);
end


for xxx = 1:length(fractalIDset)
    Regression_x = [];
    All_FR = [];
    for iii = 1:length(Neuronlist_learning)
        if strcmpi(Neuronlist_learning(iii).monkeyName, 'Lemmy')
            fractalIDset = fractalIDset_Lemmy;
        elseif strcmpi(Neuronlist_learning(iii).monkeyName, 'Slayer')
            fractalIDset = fractalIDset_Slayer;
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


scatter(Regression_L1_x,All_FR_L1,5,'b', 'filled');

xlabel('Appearance times')
ylabel('Z-scored firing rate');


fitoptions = optimset('MaxIter',10000,'MaxFunEvals',20000);
exp_paras = [0,0,0];
LB = -Inf(size(exp_paras));
UB = Inf(size(exp_paras));
LB(3) = 0;
UB(3) = 10;

for xxx = 3:length(fractalIDset) % regression on the learning trial
    eval(['Regression_x = ' Regression_x_names{xxx} ';']);
    eval(['All_FR = ' All_FR_names{xxx} ';']);
    
    fun_L = @(x) exp_fun_fitting(Regression_x, x, All_FR);
    x0 = exp_paras;
    exp_paras = fminsearchbnd(fun_L,x0,LB,UB,fitoptions);%%%
    
    eval([exp_paras_names{xxx} ' = exp_paras;']);
end


AllmeanFR_NL = mean(All_FR_NL);
AllmeanFR_FL = mean(All_FR_FL);

plot_x = 1:max(Regression_L1_x);
exponentfun = @(C,t) C(1)+C(2)*exp(-C(3)*t);
plot(plot_x,repmat(AllmeanFR_NL,1,length(plot_x)),'r');hold on;
plot(plot_x,repmat(AllmeanFR_FL,1,length(plot_x)),'k');
plot(plot_x,exponentfun(exp_paras_L1,plot_x),'b');
plot(plot_x,exponentfun(exp_paras_L2,plot_x),'color',[0,0.3,0.7]);
try
    plot(plot_x,exponentfun(exp_paras_L3,plot_x),'color',[0,0.5,0.5]);
    plot(plot_x,exponentfun(exp_paras_L4,plot_x),'color',[0,0.7,0.3]);
    plot(plot_x,exponentfun(exp_paras_L5,plot_x),'color',[0,0.9,0.1]);
catch
end

xlim([0,20]);
ylim([-0.5,1]);
title(['lambta 1day=' mat2str(1/exp_paras_L1(3),3) ', lambta 2 day=' mat2str(1/exp_paras_L2(3),3)]) ;

nsubplot(75,75,1:20,25+(1:20));
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
legend('novel', 'familiar', 'learning 1 day', 'learning n>=2 day');
xlabel('Appearance times')
ylabel('Z-scored firing rate');
title([Select_criteria{xy} ', n=' mat2str(num)]) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Consecutive plots (The derivative of Firing rate changes by appearance
nsubplot(75,75,25+(1:20),(1:20));
% ROCs between novel and learning 1 day,
%learning 1 day and learning n>=2 days
% learningn>=2 days and familiar

der_ROC_FL = zeros(1,max(Regression_L1_x)-1);
der_ROC_NL = zeros(1,max(Regression_L1_x)-1);
der_ROC_L1 = zeros(1,max(Regression_L1_x)-1);
der_ROC_L2 = zeros(1,max(Regression_L1_x)-1);

for ii = 1:(max(Regression_L1_x)-1)
    der_ROC_FL(ii) = -(rocarea3(All_FR_FL(Regression_FL_x==ii), All_FR_FL(Regression_FL_x==ii+1))-0.5);
    der_ROC_NL(ii) = -(rocarea3(All_FR_NL(Regression_NL_x==ii), All_FR_NL(Regression_NL_x==ii+1))-0.5);
    der_ROC_L1(ii) = -(rocarea3(All_FR_L1(Regression_L1_x==ii), All_FR_L1(Regression_L1_x==ii+1))-0.5);
    der_ROC_L2(ii) = -(rocarea3(All_FR_L2(Regression_L2_x==ii), All_FR_L2(Regression_L2_x==ii+1))-0.5);
end

%smoothing
der_ROC_NL = smooth(der_ROC_NL,3,'moving');
der_ROC_FL = smooth(der_ROC_FL,3,'moving');
der_ROC_L1 = smooth(der_ROC_L1,3,'moving');
der_ROC_L2 = smooth(der_ROC_L2,3,'moving');

plot_x = 1:max(Regression_L1_x)-1;
hold on;
plot(plot_x, der_ROC_NL,'r');
plot(plot_x, der_ROC_FL,'k');
plot(plot_x, der_ROC_L1,'b');
plot(plot_x, der_ROC_L2,'color',[0,0.3,0.7]);
ylabel('Area under ROC')
%legend({'Familiar vs Learning 2 days', 'Learning 2 days vs 1 day', '1 day vs Novel'});
xlim([0,20]);
ylim([-0.05,0.1]);

%% Another way to measure learning speed
normalized_meanFR_L1 = (meanFR_L1-meanFR_FL)./(meanFR_NL-meanFR_FL);
normalized_meanFR_L2 = (meanFR_L2-meanFR_FL)./(meanFR_NL-meanFR_FL);

learning_rate_L1 = -(normalized_meanFR_L1(2:end)-normalized_meanFR_L1(1:end-1))./normalized_meanFR_L1(1:end-1);
learning_rate_L2 = -(normalized_meanFR_L2(2:end)-normalized_meanFR_L2(1:end-1))./normalized_meanFR_L2(1:end-1);

smoothfun = @(x) smooth(x,3,'moving');
learning_rate_L1_s = smoothfun(learning_rate_L1)';%smooth(learning_rate_L1,3,'moving');
learning_rate_L2_s = smoothfun(learning_rate_L2)';%smooth(learning_rate_L2,3,'moving');

if xy>1 || ~run_longtimebootstrapping
    nsubplot(75,75,25+(1:20),25+(1:20));
    hold on;
    plot_x = 1:max(Regression_L1_x);
    plot(plot_x, normalized_meanFR_L1,'b');
    plot(plot_x, normalized_meanFR_L2,'color',[0,0.3,0.7]);
    ylabel('normalized Firing rate');
    xlim([0,20]);
    ylim([0,1]);
    
    
    nsubplot(75,75,25+(1:20),50+(1:20));
    plot_x = 1:max(Regression_L1_x)-1;
    
    hold on;
    smoothfun = @(x) smooth(x,3,'moving');
    % plot(plot_x, learning_rate_L1,'b');
    % plot(plot_x, learning_rate_L2,'color',[0,0.3,0.7]);
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
    tic
    %% try to batch the bootstrapping on server
    
    pars.shuffling_num = shuffling_num;
    pars.fractalIDset_Slayer = fractalIDset_Slayer;
    pars.fractalIDset_Lemmy = fractalIDset_Lemmy;
    pars.maxappearance = max(Regression_L1_x);
    pars.smoothfun = smoothfun;
    
    tic
    datastruct = Learning_acrossday_bootstrapping(Neuronlist_learning,pars);
    t = toc
    
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
    
    
    nsubplot(75,75,25+(1:20),25+(1:20));
    hold on;
    plot_x = 1:maxappearance;
    shadedErrorBar(plot_x,normalized_meanFR_L1,bound_normalized_meanFR_L1,{'-b', 'LineWidth', 0.5}, 0);
    shadedErrorBar(plot_x,normalized_meanFR_L2,bound_normalized_meanFR_L2,{'color',[0,0.3,0.7], 'LineWidth', 0.5}, 0);
    ylabel('normalized Firing rate');
    xlim([0,20]);
    ylim([0,1]);
    
    
    nsubplot(75,75,25+(1:20),50+(1:20));
    plot_x = 1:maxappearance-1;
    hold on;
    shadedErrorBar(plot_x,learning_rate_L1_s,bound_learning_rate_L1_s,{'-b', 'LineWidth', 0.5}, 0);
    shadedErrorBar(plot_x,learning_rate_L2_s,bound_learning_rate_L2_s,{'color',[0,0.3,0.7], 'LineWidth', 0.5}, 0);
    ylabel('Learning rate (alpha)')
    %legend({'Familiar vs Learning 2 days', 'Learning 2 days vs 1 day', '1 day vs Novel'});
    xlim([0,20]);
    %%
    toc
end
set(gcf,'Position',[100 100 1100 800],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',[plotpath '/Learning_n_day_fitting' Select_criteria{xy} '.pdf']);
end


function sq_err = exp_fun_fitting(Regression_x, exp_paras, firingrate)
%%% calculate squared error given a set of parameters.
% Regression_x here should a column vector of x data, firing rate should be
% a column of y data, exp_paras is the current parameter(length: 3), and the
% function returns the error of the fitting function with the actual
% function.
sq_err = 0;
for iii = 1:size(Regression_x, 1)
    appearingtime = Regression_x(iii,1);
    yvalue = exp_paras(1)+exp_paras(2)*exp(-exp_paras(3)*appearingtime);
    sq_err = sq_err+(firingrate(iii)-yvalue)^2;
end
    
end