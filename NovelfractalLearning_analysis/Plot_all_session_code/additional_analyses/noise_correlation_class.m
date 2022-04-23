classdef noise_correlation_class < handle
    properties
        corrmean_all=[];
        condition_names = [];
        res_var_curves_n = [];
        res_var_curves_f = [];
        res_var_curves_1 = [];
        included_neuron = [];
        monkey_name = {};
        
        var_nov_ax = [];
        var_rand_ax = [];
        
        neuron_mean = [];
        neuron_var = [];
    end
    
    methods
        function obj = noise_correlation_class(included_neuron, init_struct)
            obj.included_neuron=included_neuron;
            if exist('init_struct', 'var')
                propertyname = fieldnames(init_struct);
                for ii = 1:numel(propertyname)
                    if isprop(obj, propertyname{ii})
                        obj.(propertyname{ii}) = init_struct.(propertyname{ii});
                    end
                end
            end
        end
        
        function cormean = neuron_correlation(obj, FR)
            % FR is a m*n matrix, where m is the number of fractals and n
            % is the number of neurons
            [cormatrix,p] = corr(FR, 'type', 'Pearson');
            %calculate the mean of the correlations of pairs of neurons
            corsum = 0;
            n = 0;
            for ii = 1:size(cormatrix,1)
                for jj = 1:ii-1
                    if ~isnan(cormatrix(ii,jj))
                        corsum = corsum+cormatrix(ii,jj);
                        n=n+1;
                    end
                end
            end
            cormean = corsum/n;
        end
        
        function variance = FR_var(obj, FR, project_ax)
            % FR is a m*n matrix, where m is the number of fractals and n
            % is the number of neurons
            % project_ax: 1*n vector
            BLmat= mean(FR,1);
            FR_new = FR-BLmat;
            variance = var(FR_new*project_ax',1);
            
        end
        
        
        function noise_correlation_singlesession(obj, Generaltask, Neuronlist)
            obj.monkey_name{end+1,1} = Neuronlist(1).monkeyName;
            Channelnames = {Neuronlist.name};
            % only include SPK
            Neuronlist_sub = Neuronlist(contains(Channelnames, 'SPK') & [Neuronlist.exclude]==0);
            if isempty(Neuronlist_sub)
                return;
            end
            
            fracsplitlogical = ones(1, size(Generaltask.Fractals,2));
            
            frac_sets = make_frac_sets(Generaltask, fracsplitlogical);
            
            % neuron types and fractal types
            neurontypes = {'noveltyneurons', 'otherneurons'};
            condition = {'novelfractal', 'familiarfractal'};
            
            if strcmpi(obj.included_neuron, 'noveltyresponsive')
                noveltyneuronind = find([Neuronlist_sub.P_pred_nov_vs_fam]<0.01);
            elseif strcmpi(obj.included_neuron, 'noveltyexcited')
                noveltyneuronind = find([Neuronlist_sub.P_pred_nov_vs_fam]<0.01 & [Neuronlist_sub.pred_nov_vs_fam]>0);
            elseif strcmpi(obj.included_neuron, 'noveltyinhibited')
                noveltyneuronind = find([Neuronlist_sub.P_pred_nov_vs_fam]<0.01 & [Neuronlist_sub.pred_nov_vs_fam]<0);
            else
                error('No such included_neuron input');
            end
            otherneuronind = find([Neuronlist_sub.P_pred_nov_vs_fam]>0.01);
            
            frac_ind_set = [];
            frac_ind_set.frac_novel_ind = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [7999]) & Generaltask.Fractals(2,:)==6));
            frac_ind_set.frac_6401_ind = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6401])));
            frac_ind_set.frac_6404_ind = intersect(frac_sets.successfulfrac_noviol, find(ismember(Generaltask.Fractals(1,:), [6404])));
            
            frac_ind_set.frac_6500_ind = intersect(frac_sets.successfulfrac_noviol, find(Generaltask.Fractals(1,:)==6500 & Generaltask.Fractals(3,:)==2));
            frac_ind_set.frac_6501_ind = intersect(frac_sets.successfulfrac_noviol, find(Generaltask.Fractals(1,:)==6501 & Generaltask.Fractals(3,:)==2));
            frac_ind_set.frac_6502_ind = intersect(frac_sets.successfulfrac_noviol, find(Generaltask.Fractals(1,:)==6502 & Generaltask.Fractals(3,:)==2));
            frac_ind_set.frac_6503_ind = intersect(frac_sets.successfulfrac_noviol, find(Generaltask.Fractals(1,:)==6503 & Generaltask.Fractals(3,:)==2));
            frac_ind_set.frac_6504_ind = intersect(frac_sets.successfulfrac_noviol, find(Generaltask.Fractals(1,:)==6504 & Generaltask.Fractals(3,:)==2));
            frac_ind_set.frac_6505_ind = intersect(frac_sets.successfulfrac_noviol, find(Generaltask.Fractals(1,:)==6505 & Generaltask.Fractals(3,:)==2));
            
            neurontypes_ind = {noveltyneuronind, otherneuronind};
            condition_ind = {frac_ind_set.frac_novel_ind, [frac_ind_set.frac_6500_ind, frac_ind_set.frac_6501_ind, frac_ind_set.frac_6502_ind, ...
                frac_ind_set.frac_6503_ind, frac_ind_set.frac_6504_ind, frac_ind_set.frac_6505_ind] };
            
            FR_all = horzcat(Neuronlist_sub.All_Fractal_FR);
            % subtracted the signal from the firing rate.
            fnames = fieldnames(frac_ind_set);
            for ii = 1:numel(fnames)
                tempind = frac_ind_set.(fnames{ii});
                FR_all(tempind,:) = FR_all(tempind,:) - mean(FR_all(tempind,:),1);
            end
            
            
            if numel(noveltyneuronind)>=5
                if strcmpi(obj.included_neuron, 'noveltyresponsive')
                    % flip the sign of negative novelty responsive neuron
                    neg_ind = find([Neuronlist_sub.pred_nov_vs_fam]<0);
                    FR_all(:, neg_ind) = -FR_all(:, neg_ind);
                end
                
                for typeind = 1:numel(neurontypes)
                    for condind = 1:numel(condition)
                        FR = FR_all(condition_ind{condind},neurontypes_ind{typeind});
                        cormean(1,(typeind-1)*numel(condition)+condind) = obj.neuron_correlation(FR);
                        % record the names
                        obj.condition_names{end+1} = strcat(neurontypes{typeind}, ' on', ' ' ,condition{condind});
                    end
                end
                obj.corrmean_all(end+1,:) = cormean;
                
                %% variance analysis.
                bootstrap_num = 10000;
                FR_all = horzcat(Neuronlist_sub(noveltyneuronind).All_Fractal_FR);
                meannov_response = mean(FR_all(frac_ind_set.frac_novel_ind,:),1);
                meanfam_response = mean(FR_all([frac_ind_set.frac_6401_ind, frac_ind_set.frac_6404_ind],:),1);
                nov_axis = (meannov_response-meanfam_response)/norm(meannov_response-meanfam_response);
                
                temp_var_nov_ax = zeros(bootstrap_num, 3);
                temp_var_rand_ax = zeros(bootstrap_num, 3);
                temp_FR_5_nov_mean = zeros(bootstrap_num, numel(noveltyneuronind));
                temp_FR_5_fam_mean = zeros(bootstrap_num, numel(noveltyneuronind));
                temp_FR_1_fam_mean = zeros(bootstrap_num, numel(noveltyneuronind));
                temp_FR_5_nov_var = zeros(bootstrap_num, numel(noveltyneuronind));
                temp_FR_5_fam_var = zeros(bootstrap_num, numel(noveltyneuronind));
                temp_FR_1_fam_var = zeros(bootstrap_num, numel(noveltyneuronind));
                for ii = 1: bootstrap_num
                    FR_n = FR_all(frac_ind_set.frac_novel_ind(randperm(numel(frac_ind_set.frac_novel_ind),5)),:);
                    FR_1 = FR_all(frac_ind_set.frac_6500_ind(randperm(numel(frac_ind_set.frac_6500_ind), 5)),:);
                    
                    uniqe_frac_ind = [];
                    for fractID = 6501:6505
                        temp_ind = intersect(frac_sets.successfulfrac_noviol, find(Generaltask.Fractals(1,:)==fractID & Generaltask.Fractals(3,:)==2));
                        uniqe_frac_ind(end+1) = temp_ind(randi(numel(temp_ind),1));
                    end
                    FR_f = FR_all(uniqe_frac_ind(1:5),:);
                    
                    % random axis as control
                    rand_axis = randn(size(nov_axis));
                    rand_axis = rand_axis/norm(rand_axis);
                    
                    temp_var_nov_ax(ii,:) = [obj.FR_var(FR_n, nov_axis), obj.FR_var(FR_f, nov_axis), obj.FR_var(FR_1, nov_axis)];
                    temp_var_rand_ax(ii,:) = [obj.FR_var(FR_n, rand_axis), obj.FR_var(FR_f, rand_axis), obj.FR_var(FR_1, rand_axis)];
                    
                    % get the mean and var for each neuron in the three
                    % cases
                    temp_FR_5_nov_mean(ii,:) = mean(FR_n,1);
                    temp_FR_5_fam_mean(ii,:) = mean(FR_f,1);
                    temp_FR_1_fam_mean(ii,:) = mean(FR_1,1);
                    temp_FR_5_nov_var(ii,:) = var(FR_n,1,1);
                    temp_FR_5_fam_var(ii,:) = var(FR_n,1,1);
                    temp_FR_1_fam_var(ii,:) = var(FR_n,1,1);
                    
                end
                obj.var_nov_ax(end+1,:) = mean(temp_var_nov_ax,1);
                obj.var_rand_ax(end+1,:) = mean(temp_var_rand_ax,1);
                
                obj.neuron_mean = vertcat(obj.neuron_mean, [mean(temp_FR_5_nov_mean,1)', mean(temp_FR_5_fam_mean,1)', mean(temp_FR_1_fam_mean,1)']);
                obj.neuron_var = vertcat(obj.neuron_var, [mean(temp_FR_5_nov_var,1)', mean(temp_FR_5_fam_var,1)', mean(temp_FR_1_fam_var,1)']);
                
                
            else
                obj.corrmean_all(end+1,:) = nan(1,4);
                obj.var_nov_ax(end+1,:) = nan(1,3);
                obj.var_rand_ax(end+1,:) = nan(1,3);
            end
        end
        
        
        function noise_correlation_analysis(obj, path, varargin)
            normal = 0; % normalizing?
            Monkey = 'Combine';
            
            tempInd = find(strcmpi(varargin, 'normalize'),1);
            if ~isempty(tempInd)
                normal = varargin{tempInd+1};
            end
            
            tempInd = find(strcmpi(varargin, 'Monkey'),1);
            if ~isempty(tempInd)
                Monkey = varargin{tempInd+1};
            end
            
            if strcmpi(Monkey, 'Combine')
                included = true(size(obj.monkey_name));%true(size(obj.corrmean_all,1),1);
            else
                included = strcmpi(obj.monkey_name, Monkey);
            end
            
            
            figure
            nsubplot(169,169, 1:70, 1:70);
            Datamean = nanmean(obj.corrmean_all(included,:),1);
            Datastd = nanstd(obj.corrmean_all(included,:),[],1)/sqrt(sum(~isnan(obj.corrmean_all(included,1))));
            bar([1,2,3,4],Datamean, 'b'); hold on;
            errorbar([1,2,3,4],Datamean, Datastd,'r.');
            
            p = signrank(obj.corrmean_all(included,1), obj.corrmean_all(included,2));
            text(1, Datamean(1)*1.2, ['p(1vs2) = ' mat2str(p,4)]);
            p = signrank(obj.corrmean_all(included,1), obj.corrmean_all(included,3));
            text(3, Datamean(3)*1.2, ['p(1vs3) = ' mat2str(p,4)]);
            p = signrank(obj.corrmean_all(included,2), obj.corrmean_all(included,4));
            text(4, Datamean(4)*1.2, ['p(2vs4) = ' mat2str(p,4)]);
            
            % combine the novelty neuron's corr and other neuron's corr and
            % compare the two.
            p = signrank(mean(obj.corrmean_all(included,1:2),2), mean(obj.corrmean_all(included,3:4),2));
            text(2, Datamean(2)*1.2, ['p(1&2 vs 3&4) = ' mat2str(p,4)]);
            
            set(gca, 'xtick', [1,2,3,4],'xticklabel', obj.condition_names, 'xticklabelRotation', -45)
            title(['Noise correlation, n = ' mat2str(sum(~isnan(obj.corrmean_all(included,1))))]);
            
            
            set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
            print(gcf,'-dpdf', '-painters',fullfile(path,['noise_correlation_analysis_' obj.included_neuron '.pdf']));
            
        end
        
        
        function noise_variance_analysis(obj, path)
            %% variance analysis
            figure;
            nsubplot(169,169, 10:70, 10:70);
            Datamean = nanmean(obj.var_nov_ax,1);
            Datastd = nanstd(obj.var_nov_ax,[],1)/sqrt(sum(~isnan(obj.var_nov_ax(:,1))));
            bar([1,2,3],Datamean, 'b'); hold on;
            errorbar([1,2,3],Datamean, Datastd,'r.');
            title('variance along novelty axis');
            set(gca, 'xtick', [1,2,3], 'xticklabel', {'5 Novel', '5 Familiar', '1 Familiar'}, 'xticklabelRotation', -45);
            p = signrank(obj.var_nov_ax(:,1), obj.var_nov_ax(:,2));
            text(2.2, Datamean(2), mat2str(p,4));
            p = signrank(obj.var_nov_ax(:,1), obj.var_nov_ax(:,3));
            text(3.2, Datamean(3), mat2str(p,4));
            
            nsubplot(169,169, 90:150, 10:70);
            Datamean = nanmean(obj.var_rand_ax,1);
            Datastd = nanstd(obj.var_rand_ax,[],1)/sqrt(sum(~isnan(obj.var_rand_ax(:,1))));
            bar([1,2,3],Datamean, 'b'); hold on;
            errorbar([1,2,3],Datamean, Datastd,'r.');
            title('variance along random axis');
            set(gca, 'xtick', [1,2,3], 'xticklabel', {'5 Novel', '5 Familiar', '1 Familiar'}, 'xticklabelRotation', -45);
            p = signrank(obj.var_nov_ax(:,1), obj.var_nov_ax(:,2));
            text(2.2, Datamean(2), mat2str(p,4));
            p = signrank(obj.var_nov_ax(:,1), obj.var_nov_ax(:,3));
            text(3.2, Datamean(3), mat2str(p,4));
            
            nsubplot(169,169, 10:70, 90:150);
            Datamean = nanmean(obj.var_nov_ax./obj.var_rand_ax,1);
            Datastd = nanstd(obj.var_nov_ax./obj.var_rand_ax,[],1)/sqrt(sum(~isnan(obj.var_rand_ax(:,1))));
            bar([1,2,3],Datamean, 'b'); hold on;
            errorbar([1,2,3],Datamean, Datastd,'r.');
            title('difference or ratio in variance along novel axis and random axis');
            set(gca, 'xtick', [1,2,3], 'xticklabel', {'5 Novel', '5 Familiar', '1 Familiar'}, 'xticklabelRotation', -45);
            p = signrank(obj.var_nov_ax(:,1), obj.var_nov_ax(:,2));
            text(2.2, Datamean(2), mat2str(p,4));
            p = signrank(obj.var_nov_ax(:,1), obj.var_nov_ax(:,3));
            text(3.2, Datamean(3), mat2str(p,4));
            
            
            set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
            print(gcf,'-dpdf', '-painters',fullfile(path,['noise_varance_analysis' obj.included_neuron '.pdf']));
            
            
        end
        
        
    end
    
end