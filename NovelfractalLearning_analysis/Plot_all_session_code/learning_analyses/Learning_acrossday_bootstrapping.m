% bootstrapping of Learning appearance by appearance plot.
function datastruct = Learning_acrossday_bootstrapping(Neuronlist_learning,pars)
shuffling_num = pars.shuffling_num;
fractalIDset_S = pars.fractalIDset_S;
fractalIDset_L = pars.fractalIDset_L;
smoothfun = pars.smoothfun;
maxappearance = pars.maxappearance;
%% batch the bootstrapping on server


normalized_meanFR_L1_bootstrapping = zeros(shuffling_num, maxappearance);
normalized_meanFR_L2_bootstrapping = zeros(shuffling_num, maxappearance);

learning_rate_L1_s_bootstrapping = zeros(shuffling_num, maxappearance-1);
learning_rate_L2_s_bootstrapping = zeros(shuffling_num, maxappearance-1);



parfor xyz = 1:shuffling_num
    % initialize
    Regression_FL_x = [];
    Regression_NL_x = [];
    Regression_L1_x = [];
    Regression_L2_x = [];
    All_FR_FL = [];
    All_FR_NL = [];
    All_FR_L1 = [];
    All_FR_L2 = [];
    
    Randind = randi(numel(Neuronlist_learning),[numel(Neuronlist_learning),1]);
    
    for xxx = 1:4% numel(fractalIDset)
        Regression_x = [];
        All_FR = [];
        for iii = 1:length(Neuronlist_learning)
            if strcmpi(Neuronlist_learning(Randind(iii)).monkeyName, 'L')
                fractalIDset = fractalIDset_L;
            elseif strcmpi(Neuronlist_learning(Randind(iii)).monkeyName, 'S')
                fractalIDset = fractalIDset_S;
            end
            
            for ii = 1:length(fractalIDset{xxx})
                if isfield(Neuronlist_learning(Randind(iii)).learning, ['FR' mat2str(fractalIDset{xxx}(ii))])
                    Regression_x = [Regression_x; (1:length(Neuronlist_learning(Randind(iii)).learning.(['FR' mat2str(fractalIDset{xxx}(ii))])))'];
                    All_FR = [All_FR; Neuronlist_learning(Randind(iii)).learning.(['FR' mat2str(fractalIDset{xxx}(ii))])];
                end
            end
        end
        
        
        switch xxx
            case 1
                Regression_FL_x = Regression_x;
                All_FR_FL = All_FR;
            case 2
                Regression_NL_x = Regression_x;
                All_FR_NL = All_FR;
            case 3
                Regression_L1_x = Regression_x;
                All_FR_L1 = All_FR;
            case 4
                Regression_L2_x = Regression_x;
                All_FR_L2 = All_FR;
        end
    end
    
    %initialize
    meanFR_FL = zeros(1,maxappearance);
    meanFR_NL = zeros(1,maxappearance);
    meanFR_L1 = zeros(1,maxappearance);
    meanFR_L2 = zeros(1,maxappearance);
    
    for ii = 1:maxappearance
        meanFR_FL(ii) = mean(All_FR_FL(Regression_FL_x==ii));
        meanFR_NL(ii) = mean(All_FR_NL(Regression_NL_x==ii));
        meanFR_L1(ii) = mean(All_FR_L1(Regression_L1_x==ii));
        meanFR_L2(ii) = mean(All_FR_L2(Regression_L2_x==ii));
    end
    
    %
    normalized_meanFR_L1_local = (meanFR_L1-meanFR_FL)./(meanFR_NL-meanFR_FL);
    normalized_meanFR_L2_local = (meanFR_L2-meanFR_FL)./(meanFR_NL-meanFR_FL);
    
    learning_rate_L1 = -(normalized_meanFR_L1_local(2:end)-normalized_meanFR_L1_local(1:end-1))./normalized_meanFR_L1_local(1:end-1);
    learning_rate_L2 = -(normalized_meanFR_L2_local(2:end)-normalized_meanFR_L2_local(1:end-1))./normalized_meanFR_L2_local(1:end-1);
    
    normalized_meanFR_L1_bootstrapping(xyz,:) = normalized_meanFR_L1_local;
    normalized_meanFR_L2_bootstrapping(xyz,:) = normalized_meanFR_L2_local;
    learning_rate_L1_s_bootstrapping(xyz,:) = smoothfun(learning_rate_L1);
    learning_rate_L2_s_bootstrapping(xyz,:) = smoothfun(learning_rate_L2);
    
end

datastruct.normalized_meanFR_L1_bootstrapping = normalized_meanFR_L1_bootstrapping;
datastruct.normalized_meanFR_L2_bootstrapping = normalized_meanFR_L2_bootstrapping;
datastruct.learning_rate_L1_s_bootstrapping = learning_rate_L1_s_bootstrapping;
datastruct.learning_rate_L2_s_bootstrapping = learning_rate_L2_s_bootstrapping;
end
