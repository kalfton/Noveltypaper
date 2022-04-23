function [neuron_counts, segmentmap] = neuron_counts_region_V2(segmentmap, datastruct, pars)

segmentNames = segmentmap(:,1);
cells_areas = [datastruct.(pars.anatomyCase)];

%% count the total # of cells per each area
for segmentInd = 1:numel(segmentNames)
    currentCells = ismember(cells_areas, segmentmap{segmentInd,2});
    allRegionsCells(segmentInd,1) = numel(find(currentCells == 1));
    clear currentCells;
end
if pars.applyCellThreshold
    goodRegions = find(allRegionsCells >= pars.cellThreshold & ~contains(segmentNames, pars.exclude_region));
%     include_logical = false(size(segmentNames));
%     for ii = 1:numel(segmentNames)
%         include_logical(ii) = any(strcmpi(segmentNames{ii}, pars.include_region));
%     end
%     goodRegions = include_logical;
    segmentmap = segmentmap(goodRegions,:);
    segmentNames = segmentNames(goodRegions,:);
    allRegionsCells = allRegionsCells(goodRegions);
end


experimentConditions = pars.experimentConditions ;
experimentindices = pars.experimentindices;
StatisticalThreshold = pars.sigThres;
celltype = pars.celltype;
clear significants significantCounts significantPerc pouts

%% count the significant cells in the whole brain
% significants saves the index of significant cells.
for celltypeInd = 1:numel(celltype)
    for condInd = 1:numel(pars.experimentConditions) %
        conditionname = experimentConditions{condInd};
        if strcmpi(celltype{celltypeInd}, 'positive')
            significants.(celltype{celltypeInd}).(conditionname) = find([datastruct.(experimentConditions{condInd})] <= StatisticalThreshold & [datastruct.(experimentindices{condInd})] >0 ); % find the sig
        elseif strcmpi(celltype{celltypeInd}, 'negative')
            significants.(celltype{celltypeInd}).(conditionname) = find([datastruct.(experimentConditions{condInd})] <= StatisticalThreshold & [datastruct.(experimentindices{condInd})] <0 ); % find the sig
        elseif strcmpi(celltype{celltypeInd}, 'selective')
            significants.(celltype{celltypeInd}).(conditionname) = find([datastruct.(experimentConditions{condInd})] <= StatisticalThreshold);
        end
        
    end
end

%% count the significant cells in each area
for celltypeInd = 1:numel(celltype)
    for condInd = 1:numel(pars.counting_pairs) %
        tempnames = experimentConditions(pars.counting_pairs{condInd});
        conditionname=[sprintf('%s_',tempnames{1:end-1}),tempnames{end}];
        % intersect the conditions in counting pairs
        tempInd = 1:numel(datastruct);
        for ii = 1:numel(pars.counting_pairs{condInd})
            tempInd = intersect(tempInd,significants.(celltype{celltypeInd}).(tempnames{ii}));
        end
        
        tempregionInd = [datastruct(tempInd).(pars.anatomyCase)];
        
        for segmentInd = 1:numel(segmentNames) % now count the number of occurence in each area
            currentSegmentNo = segmentmap{segmentInd,2};
            p = power(StatisticalThreshold, numel(pars.counting_pairs{condInd}));
            n = allRegionsCells(segmentInd,1);
            s = numel(find(ismember(tempregionInd,currentSegmentNo)));
            significantCounts.(celltype{celltypeInd})(segmentInd).(conditionname) = s;
            significantPerc.(celltype{celltypeInd})(segmentInd).(conditionname) =  s/n;
            Sided = 'Greater';
            pouts.(celltype{celltypeInd})(segmentInd).(conditionname) = myBinomTest(s,n,p,Sided);
            
            % for debugging
            
            %         pouts.sig(segmentInd,condInd) = myBinomTest(s,n,p,Sided);
            clear temp;
        end
    end
end

neuron_counts.allRegionsCells = allRegionsCells;
neuron_counts.significantCounts = significantCounts;
neuron_counts.significantPerc = significantPerc;
neuron_counts.pouts = pouts;




end