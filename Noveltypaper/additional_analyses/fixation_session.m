classdef fixation_session < handle
    properties
        reactiontime_singletrials = cell(1,4);
        reactiontime_types = [];
        reactionpupil_types = [];
    end
    
    methods
        function obj = fixation_session(init_struct)
            if exist('init_struct', 'var')
                propertyname = fieldnames(init_struct);
                for ii = 1:numel(propertyname)
                    if isprop(obj, propertyname{ii})
                        obj.(propertyname{ii}) = init_struct.(propertyname{ii});
                    end
                end
            end
        end
        
        function fixation_data_collect(obj, Generaltask)
            % reaction time
            successtrial = [Generaltask.successtrial];
            
            Generaltask.plottype = Generaltask.trialtype;
            Generaltask.plottype(Generaltask.trialtype==6) = 1;
            Generaltask.plottype(ismember(Generaltask.trialtype,[2,3])) = 2;
            Generaltask.plottype(ismember(Generaltask.trialtype,[4,5])) = 3;
            Generaltask.plottype(Generaltask.trialtype==1) = 4;
            
            % reaction time
            reactiontimes = Generaltask.fixacq(:,1) - Generaltask.timefpon(:,1);
            loc_reactiontime_types = zeros(4,1);
            for trialtype = 1:4
                temp_times = reactiontimes(Generaltask.plottype ==trialtype & successtrial);
                loc_reactiontime_types(trialtype) = mean(temp_times);
                obj.reactiontime_singletrials{trialtype} = [obj.reactiontime_singletrials{trialtype}; temp_times];
            end
            obj.reactiontime_types(1:4,end+1) = loc_reactiontime_types;
            
        end
        
        function fixation_analysis(obj, path)
            % path: where to save the figures
            if nargin<2
                path = '';
            end
            % reaction time
            figure;
            nsubplot(2,2,1,1);
            Datamean = mean(obj.reactiontime_types,2);
            Datastd = std(obj.reactiontime_types,[],2)/sqrt(size(obj.reactiontime_types,2));
            bar([1,2,3,4],Datamean, 'b'); hold on;
            errorbar([1,2,3,4],Datamean, Datastd,'r.');

            p = signrank(obj.reactiontime_types(1,:), obj.reactiontime_types(2,:));
            text(1, Datamean(1)*1.2, ['p(1vs2) = ' mat2str(p,4)]);
            p = signrank(obj.reactiontime_types(2,:), obj.reactiontime_types(3,:));
            text(2, Datamean(2)*1.2, ['p(2vs3) = ' mat2str(p,4)]);
            p = signrank(obj.reactiontime_types(3,:), obj.reactiontime_types(1,:));
            text(3, Datamean(3)*1.2, ['p(3vs1) = ' mat2str(p,4)]);
            ylim([0.1, 0.2]);
            set(gca, 'xtick', [1,2,3,4])
            title(['Reaction time to Fixation(session by session)']);
            
            nsubplot(2,2,1,2);
            Datamean = cellfun(@(x) mean(x),obj.reactiontime_singletrials);
            Datastd = cellfun(@(x) std(x)/sqrt(numel(x)),obj.reactiontime_singletrials);
            bar([1,2,3,4],Datamean, 'b'); hold on;
            errorbar([1,2,3,4],Datamean, Datastd,'r.');
            p = ranksum(obj.reactiontime_singletrials{1}, obj.reactiontime_singletrials{2});
            text(1, Datamean(1)*1.2, ['p(1vs2) = ' mat2str(p,4)]);
            p = ranksum(obj.reactiontime_singletrials{2}, obj.reactiontime_singletrials{3});
            text(2, Datamean(2)*1.2, ['p(2vs3) = ' mat2str(p,4)]);
            p = ranksum(obj.reactiontime_singletrials{3}, obj.reactiontime_singletrials{1});
            text(3, Datamean(3)*1.2, ['p(3vs1) = ' mat2str(p,4)]);
            ylim([0.1, 0.2]);
            set(gca, 'xtick', [1,2,3,4])
            title(['Reaction time to Fixation (trial by trial)']);
            print(fullfile(path,['Reactiontime_2_fixation.pdf']), '-dpdf')
            
        end
        
        
    end
    
end