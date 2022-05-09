function [area,pval,permarea] = rocarea3(N,S,pval_method)
% [area,pval,permarea] = rocarea3(N,S,pval_method)
%
% speedier version of rocarea, works in 2 dimensions, as follows:
%  requires S(ignal) and N(oise) to be matrices with the same number of columns
%  
%  for each pair of columns S(:,c) and N(:,c), 
%   calculates ROC for discriminating S from N:
%    0.0 if all N > all S
%    0.5 if not discriminable by a fixed criterion
%    1.0 if all S > all N
%
%  if specify nperms as a scalar, does that many permutations
%   and returns the 2-tailed p-value for each of the calculated ROC areas
%   and, if requested, returns all of the permuted-data ROCs
%   in a (nperms) x (ncolumns) matrix
%
%  pval_method can be:
% - 'approximate' or []: uses the 'ranksum' function's approximate wilcoxon test
%     based on a normal approximation to its sampling distribution.
%     NOTE: for high ROCs, this will underestimate the p-value by several
%     orders of magnitude! (this is the default)
% - a scalar: use a permutation test, doing the number of
%    specified permutations. if the output argument 'permarea' is 
%    specified, also saves the ROC areas computed from permuted data, 
%    and returns them.
% - 'exact': uses the 'ranksum' function's exact test, using the network
%    algorithm of Mehta, Patel and Tsiatis (1984).
%    NOTE: this can be rather slow, especially with a moderate or greater
%    number of data points (>30). However, if there is little overlap
%    between the two distributions (e.g. area > .9), it can be very fast;
%    if the roc area is 1.0, the calculation is near instantaneous!
% - 'default': uses 'exact' for N < 20, 'approximate' for N >= 20

assert(size(N,2) == size(S,2),'N and S must have same # columns');

if nargin < 3
    pval_method = 'approximate';
end;
if isscalar(pval_method)
    nperms = pval_method;
end;

ns = size(S,1);
nn = size(N,1);
nc = size(S,2);
area = nan*ones(1,nc);

if ns == 0 || nn == 0
    area(:) = .5;
    pval = 1;
    permarea = [];
    return;
end;
if ns == 1 && nn == 1
    area = (sign(S-N)./2)+.5;
    pval = 1;
    permarea = [];
    if isscalar(pval_method)
        permarea = repmat(area,nperms,1);
    end;
    return;
end;

% exhaustively scan all possible 
%  discriminator criteria
N = sort(N,1);
S = sort(S,1);
pct_noise = ones(ns+nn+2,1);
pct_signal = ones(ns+nn+2,1);
max_i_per_col = nn - sum(isnan(N),1);
max_j_per_col = ns - sum(isnan(S),1);
noise_increment = 1./max_i_per_col;
signal_increment = 1./max_j_per_col;
% max_i_per_col = nn - zeros(1,nc);
% max_j_per_col = ns - zeros(1,nc);
for c = 1:nc
    i = 1;
    j = 1;
    pct_noise(:) = 1;
    pct_signal(:) = 1;
    pct_noise(1) = 0;
    pct_signal(1) = 0;
    cum_noise = 0;
    cum_signal = 0;
    k = 2;
    curcrit = min(N(i,c),S(j,c));
    while i <= max_i_per_col(c) && j <= max_j_per_col(c)
        if curcrit==N(i,c)
            cum_noise = cum_noise + noise_increment(c);
            i = i +1;
        elseif curcrit==S(j,c)
            cum_signal = cum_signal + signal_increment(c);
            j = j +1;
        else
            pct_noise(k) = cum_noise;
            pct_signal(k) = cum_signal;
            k = k + 1;

            curcrit = min(N(i,c),S(j,c));
        end;
    end;
    while j<=max_j_per_col(c) && curcrit==S(j,c)
        cum_signal = cum_signal + signal_increment(c);
        j = j + 1;
    end;
    while i<=max_i_per_col(c) && curcrit==N(i,c)
        cum_noise = cum_noise + noise_increment(c);
        i = i + 1;
    end;
    pct_noise(k) = cum_noise;
    pct_signal(k) = cum_signal;
    area(c) = sum(diff(pct_signal).*(pct_noise(1:(end-1)) + 0.5.*diff(pct_noise)));
end;


if nargout > 1
    pval = nan*ones(1,nc);
    switch pval_method
        case {'approximate','exact','default'}
            if isequal(pval_method,'default')
                if ns+nn < 20
                    pval_method = 'exact';
                else
                    pval_method = 'approximate';
                end;
            end;
            for c = 1:nc
                pval(c) = ranksum(N(:,c),S(:,c),'method',pval_method);
            end;
            % pval might be set to NaN by the 'approximate' method
            %  if all of the Noise and Signal values were equal to each 
            %  other. In that case, the p-value should be 1
            pval(isnan(pval)) = 1;
        otherwise
            assert(isscalar(pval_method),'pval_method must be ''approximate'' ''exact'' or a scalar');
            permarea = nan*ones(nperms,nc);
            if nargin < 3 || isempty(nperms) || nperms < 1
                return;
            end;
            A = [N ; S];
            for p = 1:nperms
                Anew = A(randperm(ns+nn),:);
                permarea(p,:) = rocarea2(Anew(1:nn,:),Anew((nn+1):end,:));
            end;
            for c = 1:nc
                left_tail_pval = sum(permarea(:,c) <= area(c)) ./ nperms;
                right_tail_pval = sum(permarea(:,c) >= area(c)) ./ nperms;
                pval(c) = min(1,2*min(left_tail_pval,right_tail_pval));
            end;
    end;
end;

