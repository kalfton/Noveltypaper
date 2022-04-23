classdef novelty_classifier_class < handle
    properties
        train_correct_rate_all = [];
        fbeta_train_all = [];
        test_correct_rate_all = [];
        fbeta_test_all = [];
        condition_names = [];
    end
    
    methods
        
        function obj = novelty_classifier_class(init_struct)
            if exist('init_struct', 'var')
                propertyname = fieldnames(init_struct);
                for ii = 1:numel(propertyname)
                    if isprop(obj, propertyname{ii})
                        obj.(propertyname{ii}) = init_struct.(propertyname{ii});
                    end
                end
            end
        end
        
        function novelty_classifier_singlesession(obj, Generaltask, Neuronlist, responsive)
            % this function goes through all files and it will train and test classifiers for each file.
            % responsive should be 0 or 1, 0 means include all neurons, 1
            % means only include xxx-responsive neurons
            if ~exist('responsive', 'var')
                responsive = 0;
            end
            Channelnames = {Neuronlist.name};
            % only include SPK
            %Channelnames = Channelnames(contains(Channelnames, 'SPK'));
            Neuronlist_sub = Neuronlist(contains(Channelnames, 'SPK')& [Neuronlist.exclude]==0);
            
            
            fracsplitlogical = ones(1, size(Generaltask.Fractals,2));
            
            
            frac_sets = make_frac_sets(Generaltask, fracsplitlogical);
            Generaltask.Fractals = IFI(Generaltask.Fractals, frac_sets);
            frac_sets = make_recency_frac_sets(Generaltask.Fractals,frac_sets);
            
            % 4 conditions {training, testing}
            conditions = {{'surprise', 'novelty'}, {'recency', 'novelty'}, {'novelty','surprise'}, {'novelty','recency'}};
            
            train_correct_rate = zeros(1, numel(conditions));
            fbeta_train = zeros(1, numel(conditions));
            test_correct_rate = zeros(1, numel(conditions));
            fbeta_test = zeros(1, numel(conditions));
            
            for condind = 1: numel(conditions)
                switch conditions{condind}{1}
                    case 'recency'
                        train_ind = [frac_sets.shortIFI_ind, frac_sets.longIFI_ind];
                        train_y = Generaltask.Fractals(7,train_ind)'<1.5; % shortITI 1, longITI 0
                        % balancing the train label
                        tempind1 = find(train_y==0);
                        tempind2 = find(train_y==1);
                        sampling_n = min([numel(tempind1), numel(tempind2)]);
                        tempind1 = randsample(tempind1,sampling_n);
                        tempind2 = randsample(tempind2,sampling_n);
                        train_ind = train_ind([tempind1;tempind2]);
                        train_y = train_y([tempind1;tempind2]);
                        
                    case 'surprise'
                        train_ind = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6400:6510]) & Generaltask.Fractals(3,:)==3) );
                        train_y = Generaltask.Fractals(2,train_ind)'<3.5; % trialtype 2,3 predictable, trialtype 4,5 unpredictable.
                        % balancing the train label
                        tempind1 = find(train_y==0);
                        tempind2 = find(train_y==1);
                        sampling_n = min([numel(tempind1), numel(tempind2)]);
                        tempind1 = randsample(tempind1,sampling_n);
                        tempind2 = randsample(tempind2,sampling_n);
                        train_ind = train_ind([tempind1;tempind2]);
                        train_y = train_y([tempind1;tempind2]);
                        
                    case 'novelty'
                        train_ind = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6400:6410]) & Generaltask.Fractals(3,:)==2 ...
                            | Generaltask.Fractals(1,:)==7999 & Generaltask.Fractals(2,:) == 6) );
                        train_y = Generaltask.Fractals(2,train_ind)'<5.5; % trialtype 2,3 predictable familiar, trialtype 6 novel.
                        % balancing the train label
                        tempind1 = find(train_y==0);
                        tempind2 = find(train_y==1);
                        sampling_n = min([numel(tempind1), numel(tempind2)]);
                        tempind1 = randsample(tempind1,sampling_n);
                        tempind2 = randsample(tempind2,sampling_n);
                        train_ind = train_ind([tempind1;tempind2]);
                        train_y = train_y([tempind1;tempind2]);
                        
                    otherwise
                        error("invalid input");
                end
                
                switch conditions{condind}{2}
                    case 'recency'
                        test_ind = [frac_sets.shortIFI_ind, frac_sets.longIFI_ind];
                        test_y = Generaltask.Fractals(7,test_ind)'<1.5; % shortITI 1, longITI 0
                        % balancing the test label
                        tempind1 = find(test_y==0);
                        tempind2 = find(test_y==1);
                        sampling_n = min([numel(tempind1), numel(tempind2)]);
                        tempind1 = randsample(tempind1,sampling_n);
                        tempind2 = randsample(tempind2,sampling_n);
                        test_ind = test_ind([tempind1;tempind2]);
                        test_y = test_y([tempind1;tempind2]);
                        
                    case 'surprise'
                        test_ind = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6400:6510]) & Generaltask.Fractals(3,:)==3) );
                        test_y = Generaltask.Fractals(2,test_ind)'<3.5; % trialtype 2,3 predictable, trialtype 4,5 unpredictable.
                        % balancing the test label
                        tempind1 = find(test_y==0);
                        tempind2 = find(test_y==1);
                        sampling_n = min([numel(tempind1), numel(tempind2)]);
                        tempind1 = randsample(tempind1,sampling_n);
                        tempind2 = randsample(tempind2,sampling_n);
                        test_ind = test_ind([tempind1;tempind2]);
                        test_y = test_y([tempind1;tempind2]);
                        
                    case 'novelty'
                        test_ind = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6400:6410]) & Generaltask.Fractals(3,:)==2 ...
                            | Generaltask.Fractals(1,:)==7999 & Generaltask.Fractals(2,:) == 6) );
                        test_y = Generaltask.Fractals(2,test_ind)'<5.5; % trialtype 2,3 predictable familiar, trialtype 6 novel.
                        % balancing the test label
                        tempind1 = find(test_y==0);
                        tempind2 = find(test_y==1);
                        sampling_n = min([numel(tempind1), numel(tempind2)]);
                        tempind1 = randsample(tempind1,sampling_n);
                        tempind2 = randsample(tempind2,sampling_n);
                        test_ind = test_ind([tempind1;tempind2]);
                        test_y = test_y([tempind1;tempind2]);
                        
                    otherwise
                        error("invalid input");
                end
                
                
                switch conditions{condind}{1}
                    case 'recency'
                        responsive_ind = [Neuronlist_sub.P_recency_ind_match_pos]<0.01;
                    case 'surprise'
                        responsive_ind = [Neuronlist_sub.P_pred_vs_unpred_fam_perm]<0.01;
                    case 'novelty'
                        responsive_ind = [Neuronlist_sub.P_pred_nov_vs_fam]<0.01;
                    otherwise
                        error("invalid input");
                end
                if responsive
                    Neuronlist_resp = Neuronlist_sub(responsive_ind);
                else
                    Neuronlist_resp = Neuronlist_sub;
                end
                
                
                if numel(Neuronlist_resp)>=5 % only train the sessions with more than 5 neurons
                    train_x = [];
                    test_x = [];
                    for ii = 1: numel(Neuronlist_resp)
                        %use SVM to classify the fractals
                        train_x(end+1,:) = Neuronlist_resp(ii).All_Fractal_FR(train_ind);
                        test_x(end+1,:) = Neuronlist_resp(ii).All_Fractal_FR(test_ind);
                    end
                    train_x = train_x';
                    test_x = test_x';
                    
                    
                    SVMModel = fitcsvm(train_x, train_y);
                    % evaluate on train set
                    [train_predict_y,score] = predict(SVMModel, train_x);
                    % correct rate
                    train_correct_rate(1,condind) = sum(train_y == train_predict_y)/numel(train_y);
                    fbeta_train(1,condind) = fbeta_score(train_predict_y, train_y);
                    
                    % evaluate on test set
                    [test_predict_y,score] = predict(SVMModel, test_x);
                    % correct rate
                    test_correct_rate(1,condind) = sum(test_y == test_predict_y)/numel(test_y);
                    fbeta_test(1,condind) = fbeta_score(test_predict_y, test_y);
                else
                    train_correct_rate(1,condind) = nan;
                    fbeta_train(1,condind) = nan;
                    test_correct_rate(1,condind) = nan;
                    fbeta_test(1,condind) = nan;
                end
                
            end
            
            obj.train_correct_rate_all(end+1,:) = train_correct_rate;
            obj.fbeta_train_all(end+1,:) = fbeta_train;
            obj.test_correct_rate_all(end+1,:) = test_correct_rate;
            obj.fbeta_test_all(end+1,:) = fbeta_test;
            obj.condition_names = cellfun(@(x) strjoin(x, '-> '), conditions, 'UniformOutput',false);
            
        end
        
        function novelty_classifier_analysis(obj, path)
            %training correct rate
            correct_thr = 0.5;
            figure;
            nsubplot(169,169, 1:70, 1:70);
            Datamean = nanmean(obj.train_correct_rate_all,1);
            Datastd = nanstd(obj.train_correct_rate_all,[],1)./sqrt(sum(~isnan(obj.train_correct_rate_all),1));
            bar([1,2,3,4],Datamean, 'b'); hold on;
            errorbar([1,2,3,4],Datamean, Datastd,'r.', 'CapSize',10);
            plot(get(gca, 'xlim'),[0.5,0.5], 'k');
            for ii = 1:4
                p = signrank(obj.train_correct_rate_all(:,ii),0.5);
                text(ii,1, ['p = ' mat2str(p,4)]);
            end
            set(gca, 'xtick', [1,2,3,4], 'xticklabel', obj.condition_names)
            ylabel('Correction rate');
            title(['classifier training correction rate']);
            
            %testing correct rate
            nsubplot(169,169, 81:150, 81:150);
            % statistical test
            Datamean = zeros(1,4);
            Datastd = zeros(1,4);
            for ii = 1:4
                include_ind = obj.train_correct_rate_all(:,ii)>correct_thr;
                Datamean(ii) = mean(obj.test_correct_rate_all(include_ind,ii));
                Datastd(ii) = std(obj.test_correct_rate_all(include_ind,ii))./sqrt(sum(include_ind));
                p = signrank(obj.test_correct_rate_all(include_ind,ii)-0.5);
                text(ii,0.65, ['p = ' mat2str(p,4)]);
            end
            
            bar([1,2,3,4],Datamean, 'b'); hold on;
            errorbar([1,2,3,4],Datamean, Datastd,'r.', 'CapSize',10);
            plot(get(gca, 'xlim'),[0.5,0.5], 'k');
            set(gca, 'xtick', [1,2,3,4], 'xticklabel', obj.condition_names)
            
            ylim([0.4,0.65]);
            ylabel('Correction rate');
            title(['classifier testing correction rate']);
            set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
            print(gcf,'-dpdf', '-painters',fullfile(path,['novelty_classifier_analysis.pdf']));
        end
        
        
    end
end

