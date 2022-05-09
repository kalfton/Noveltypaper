function [pval,r] = permutation_pair_test_fast(x,y,nperms,type)
% [pval,r] = permutation_pair_test_fast(x,y,nperms,type)
%
% calculate p-value of the test statistic: r(x,y)
%  using the null hypothesis that the x's and y's are drawn independently
%  of each other.
%
% custom written fast permutation test, for specific test functions:
%
%  'rankcorr' - Spearman rank correlation, test for difference from 0
%  'corr' - Pearson linear correlation, test for difference from 0
%  'regress' - regression coefficient b(2) in equation y = b(1) + b(2)*x + N(0,sigma),
%              test for difference from 0
%

% remove NaNs
x = x(:);
y = y(:);
ok = ~(isnan(x) | isnan(y));
x = x(ok);
y = y(ok);

n = numel(x);
assert(numel(x) == numel(y));

switch type
    % 'rankcorr' - SPEARMAN RANK CORRELATION
    case 'rankcorr'
        if n < 2
            r = nan;
            pval = 1;
            return;
        end;

        % calculate rank correlation
        n3const = (n+1)*n*(n-1) ./ 3;

        % real data
        [xrank, xadj] = fast_tiedrank(x,n);
        [yrank, yadj] = fast_tiedrank(y,n);
        D = sum((xrank - yrank).^2);
        meanD = (n3const - (xadj+yadj)./3) ./ 2;
        stdD = sqrt((n3const./2 - xadj./3)*(n3const./2 - yadj./3)./(n-1));
        r = (meanD - D) ./ (sqrt(n-1)*stdD);
        r(r > 1) = 1;
        r(r < -1) = -1;

        % permutations
        permr = zeros(nperms,1);
        for p = 1:nperms
            % permute data
            [ignore,pid] = sort(rand(1,n));
            yrank = yrank(pid);
            
            % recalculate correlation
            D = sum((xrank - yrank).^2);
            permr(p) = (meanD - D) ./ (sqrt(n-1)*stdD);
        end;
        permr(permr > 1) = 1;
        permr(permr < -1) = -1;


        % calculate p-value
        pval = min(1, 2*min(mean(permr <= r), mean(permr >= r)) );
        
    % 'corr' - PEARSON LINEAR CORRELATION
    case 'corr'
        if n < 2
            r = nan;
            pval = 1;
            return;
        end;
        
        % corr(x,y) = cov(x,y) ./ sqrt(var(x)*var(y))
        
        meanx = mean(x);
        meany = mean(y);
        xres = x - meanx;
        yres = y - meany;

        varx = (xres' * xres) ./ (n-1);
        vary = (yres' * yres) ./ (n-1);
        covxy = (xres' * yres) ./ (n-1);
        r = covxy ./ sqrt(varx .* vary);
        
        permr = zeros(nperms,1);
        for p = 1:nperms
            % shuffle xresiduals and recalculate (xres * yres) term.
            [ignore,pid] = sort(rand(1,n));
            xres = xres(pid);
            permr(p) = (xres' * yres);
        end;
        % convert to correlation
        permr = permr ./ (n-1);
        permr = permr ./ sqrt(varx .* vary);
        
        % calculate p-value
        pval = min(1, 2*min(mean(permr <= r), mean(permr >= r)) );
        
    case 'regress'
        if n < 2
            r = nan;
            pval = 1;
            return;
        end;
        
        % p-value is same as for Pearson linear correlation
        pval = permutation_pair_test_fast(x,y,nperms,'corr');
        % return regression coefficient
        r = regress(y(:),[ones(size(x(:))) x(:)]);
        r = r(2);
        
    otherwise
        error('unknown test type');
end;


% ----
function [r,tieadj] = fast_tiedrank(x,n)
tieadj = 0;

[sx, rowidx] = sort(x); 
ranks = 1:n;

ties = diff(sx)==0;
tieloc = [find(ties); n+2];
maxTies = numel(tieloc);

tiecount = 1;
while (tiecount < maxTies)
    tiestart = tieloc(tiecount);
    ntied = 2;
    while(tieloc(tiecount+1) == tieloc(tiecount)+1)
        tiecount = tiecount+1;
        ntied = ntied+1;
    end

    tieadj = tieadj + ntied*(ntied-1)*(ntied+1)/2;
    
    % Compute mean of tied ranks
    ranks(tiestart:tiestart+ntied-1) = ...
                  sum(ranks(tiestart:tiestart+ntied-1)) / ntied;
    tiecount = tiecount + 1;
end

% Broadcast the ranks back out, including NaN where required.
r(rowidx) = ranks;
