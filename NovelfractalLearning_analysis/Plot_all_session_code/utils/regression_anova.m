%% this function is aimed to use n-way anova to regress the position, trial type, fractal id etc. etc. effects.
function [FR_subtracted, anova_result]= regression_anova(Firing_rates, Fractal_info, Factor_dim, subtract_factor)
    % Factor_dim: a cell array containing the conbinations of Position,
    % fractal ID and Recency 
    % subtract_factor: array same dim as Factor_dim include the factor then
    % the input is 1, otherwise 0
    % 
    if nargin<4
        subtract_factor = ones(size(Factor_dim));
    end
    Factor_x = zeros(numel(Firing_rates), numel(Factor_dim));
    for ii = 1: numel(Factor_dim)
        factor_name = Factor_dim{ii};
        switch factor_name
            case 'position'
                Factor_x(:,ii) = Fractal_info(3,:)';
            case 'fractalID'
                Factor_x(:,ii) = Fractal_info(1,:)';
            case 'recency'
                Factor_x(:,ii) = Fractal_info(7,:)';
            otherwise
        end
    end
%     Factor_x = zeros(numel(Firing_rates), 3);
%     Factor_x(:,1) = Fractal_info(3,:)'; % order factor
%     Factor_x(:,2) = Fractal_info(1,:)'; % fractal ID factor
%     Factor_x(:,3) = Fractal_info(7,:)'; % fractal recency category
    
    
    FR_subtracted = Firing_rates; % initialize FR_subtracted
    Factor_x_origin = Factor_x;
    goodind = ~isnan(Firing_rates);% get rid of Firing rate with nan.
    Firing_rates = Firing_rates(goodind);
    Factor_x = Factor_x(goodind,:);
    
    [anova_result.anova_p,anova_result.anova_table,anova_result.anova_stats] = anovan(Firing_rates, Factor_x,'display','off','model','linear','varnames',Factor_dim);
    
    % read the factor effect and save it in eff struct
    
    for ii = find(subtract_factor)
        eff.(Factor_dim{ii}) = anova_result.anova_stats.coeffs(contains(anova_result.anova_stats.coeffnames, Factor_dim(ii)));
        % subtract the effect from raw firing rate
        
        % build a hashmap
        keySet = unique(Factor_x(:,ii));
        %valueSet = 1:length(keySet);
        subtr_val_set = eff.(Factor_dim{ii});
        if numel(keySet)~=numel(subtr_val_set)
            error('The key set and subtracted value set are different.')
        end
        %factormap = containers.Map(keySet,valueSet);
        %inds = cell2mat(values(factormap, num2cell(Factor_x(:,ii))));
        % alternative way:
        %inds = nan(size(FR_subtracted));
        subtr_val = nan(size(FR_subtracted));
        for jj = 1:numel(keySet)
            %inds(Factor_x(:,ii) == keySet(jj))=valueSet(jj);
            subtr_val(Factor_x_origin(:,ii) == keySet(jj))=subtr_val_set(jj);
        end
        %FR_subtracted_old = FR_subtracted - eff.(Factor_dim{ii})(inds);
        FR_subtracted = FR_subtracted - subtr_val;
        
    end
    

end
    