classdef sanity_check_v2 < handle
    % thoughts: define the datastruct from all sessions inside the class and
    % in each iteration of loading new session we update the datastruct
    properties
        datastruct=[]
        fieldnames = []
        Statthresh = 0.01;
    end
    
    methods
        function obj=sanity_check(inputstruct)
            if nargin>0
                obj.datastruct = inputstruct;
            else
                obj.datastruct = [];
            end
            
        end
        
        function data_collection(obj, Neuronlist)
            
            Channelnames = {Neuronlist.name};
            includedind = contains(Channelnames, 'SPK') & [Neuronlist.exclude]==0;
            
            obj.fieldnames = {'pred_vs_unpred_fam', 'P_pred_vs_unpred_fam_perm', 'pred_vs_unpred_fam_control', 'P_pred_vs_unpred_fam_perm_control', ...
                'pred_vs_unpred_fam_old', 'P_pred_vs_unpred_fam_perm_old', 'pred_nov_vs_fam_control1', 'P_pred_nov_vs_fam_control1', ...
                'pred_nov_vs_fam_control2', 'P_pred_nov_vs_fam_control2', 'pred_nov_vs_fam', 'P_pred_nov_vs_fam',...
                'obj_select_ind', 'P_object_select_anova'};
            
            
            for ii = 1: numel(obj.fieldnames)
                data(:,ii) = vertcat(Neuronlist(includedind).(obj.fieldnames{ii}));
            end
            
            
            
            
            % add data to existing datastruct
            obj.datastruct = vertcat(obj.datastruct, data);
        end
        
        function sanity_check_analysis(obj, path)
            % path: where to save the figures
            if nargin<2
                path = '';
            end
            
            plotplacesetx = {51:85,91:125,131:165,51:85,91:125,131:165,51:85,91:125,131:165, 11:45, 11:45, 11:45};
            plotplacesety = {1:49, 1:49, 1:49, 61:109, 61:109, 61:109, 121:169,121:169,121:169, 61:109, 121:169, 1:49};
            compare_xy = {[1,3], [2,4], [1,5], [2,6], [7,9], [8,10], [7,11]};
            figure;
            for xyw = 1:numel(compare_xy)
                nsubplot(169,169, plotplacesety{xyw}, plotplacesetx{xyw});
                x_data = obj.datastruct(:,compare_xy{xyw}(1));
                y_data = obj.datastruct(:,compare_xy{xyw}(2));
                nan_ind = isnan(x_data) | isnan(y_data);
                x_data(nan_ind) = [];
                y_data(nan_ind) = [];
                scatter(x_data, y_data);
                xlabel(obj.fieldnames{compare_xy{xyw}(1)}, 'Interpreter', 'none');
                ylabel(obj.fieldnames{compare_xy{xyw}(2)}, 'Interpreter', 'none');
                [rho,p] = corr(x_data, y_data, 'Type', 'Spearman');
                title(sprintf("rho = %.4f, p = %.4f", rho, p))
            end
            xyw = xyw+1;
            nsubplot(169,169, plotplacesety{xyw}, plotplacesetx{xyw});
            preselectind = obj.datastruct(:, 14)<obj.Statthresh;
            x_data = obj.datastruct(preselectind,1);
            y_data = obj.datastruct(preselectind,3);
            
            scatter(x_data, y_data);
            xlabel(obj.fieldnames{1});
            ylabel(obj.fieldnames{3});
            [rho,p] = corr(x_data, y_data, 'Type', 'Spearman');
            title(sprintf("preselected neuron by obj selectivity, rho = %.4f, p = %.4f", rho, p))
            
            set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
            print(gcf,'-dpdf', '-painters',fullfile(path,['sanity_check.pdf']));
            
        end
        
    end
    
end