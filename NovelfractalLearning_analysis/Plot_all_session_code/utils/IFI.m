%% this script is for calculating the inter-fractal interval for recency.
function Fractal_info = IFI(Fractal_info, frac_sets)
    % return Fractal info with additional rows (6): IFI, and (7): recency
    % category
    
    uniqID = unique(Fractal_info(1,:));
    uniqID = uniqID(~isnan(uniqID));
    uniqID = setdiff(uniqID, [7999]);
    n = numel(uniqID);
    Fractal_info(6:7,:) = nan;
    fractype = cell(n,1);
    timedifference = cell(n,1);
    for i=1:n
        fractype{i} = sort(intersect(frac_sets.successfulfrac,find(Fractal_info(1,:)==uniqID(i)))); 
        fracnum = fractype{i};
        timedifference{i} = Fractal_info(4,fracnum(2:end))-Fractal_info(4,fracnum(1:end-1)); % caculate time gap between the appearance of fractal
        timedifference{i} = [inf, timedifference{i}]; % the first time the fractal appears, set the time gap to be inf
        Fractal_info(6,fractype{i}) = timedifference{i};
        recency_category = ones(size(timedifference{i}));
        recency_category(timedifference{i}>2) = 2; % time difference smaller than 2s means that the previous fractal is in the same trial.
        Fractal_info(7,fractype{i}) = recency_category;
    end
    % novel fractal 7999: 
    Fractal_info(6,Fractal_info(1,:)==7999) = inf;
    Fractal_info(7,Fractal_info(1,:)==7999) = 2;
    
end