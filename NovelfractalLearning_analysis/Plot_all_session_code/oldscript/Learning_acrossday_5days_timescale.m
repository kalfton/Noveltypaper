% Fitting function within a day
tic
logical_for_neuronlist = {[Neuronlist_good(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]'>0
    [Neuronlist_good(:).P_pred_nov_vs_fam]'<StatisticalThreshold & [Neuronlist_good(:).pred_nov_vs_fam]'>0
    [Neuronlist_good(:).P_pred_nov_vs_fam]'<StatisticalThreshold};

Select_criteria = {'Novelty excited', 'Novelty inhibited', 'Novelty Other'};

for xy = 1:3 % novelty excited, inhibited and other
    

Neuronlist_learning = Neuronlist_good(logical_for_neuronlist{xy});
num = sum(logical_for_neuronlist{xy});
figure;
nsubplot(50,50,1:20,1:20);
hold on;
% initialize
for xxx = 1:length(fractalIDset)
    eval([Regression_x_names{xxx} ' = [];']);
    eval([All_FR_names{xxx} ' = [];']);
end


for xxx = 1:length(fractalIDset)
    Regression_x = [];
    All_FR = [];
    for iii = 1:length(Neuronlist_learning)
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
%scatter(Regression_FL_x,All_FR_FL,5,'k', 'filled');
%scatter(Regression_NL_x,All_FR_NL,5,'r', 'filled');
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

% fun_L = @(x) exp_fun_fitting(Regression_L_x, x, All_FR_L);
% x0 = exp_paras;
% exp_paras_L = fminsearchbnd(fun_L,x0,LB,UB,fitoptions);%%%
% 
% fun_Ln = @(x) exp_fun_fitting(Regression_Ln_x, x, All_FR_Ln);
% x0 = exp_paras;
% exp_paras_Ln = fminsearchbnd(fun_Ln,x0,LB,UB,fitoptions);%%%

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
ylim([-1,1]);
title(['lambta 1day=' mat2str(1/exp_paras_L1(3),3) ', lambta 2 day=' mat2str(1/exp_paras_L2(3),3)]) ;

nsubplot(50,50,1:20,25+(1:20));
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

plot_x = 1:max(Regression_L1_x);
hold on;
plot(plot_x, meanFR_NL,'r');
plot(plot_x, meanFR_FL,'k');
plot(plot_x, meanFR_L1,'b');
plot(plot_x, meanFR_L2,'color',[0,0.3,0.7]);
try
    plot(plot_x, meanFR_L3,'color',[0,0.5,0.5]);
    plot(plot_x, meanFR_L4,'color',[0,0.7,0.3]);
    plot(plot_x, meanFR_L5,'color',[0,0.9,0.1]);
catch
end
xlim([0,20]);
ylim([-1,1]);
legend('novel', 'familiar', 'learning 1 day', 'learning 2 day', 'learning 3 day', 'learning 4 day', 'learning 5 day');
xlabel('Appearance times')
ylabel('Z-scored firing rate');
title([Select_criteria{xy} ', n=' mat2str(num)]) ;

set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
print(gcf,'-dpdf', '-painters',[plotpath '/Learning_n_day_fitting' Select_criteria{xy} '.pdf']);
end


Time_withindayfitting = toc




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