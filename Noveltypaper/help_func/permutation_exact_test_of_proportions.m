function [pval] = permutation_exact_test_of_proportions(s1,n1,s2,n2,runmode)
% [pval] = permutation_exact_test_of_proportions(s1,n1,s2,n2,runmode)
% calculate permutation p-value of whether 
%  two proportions (s1/n1) and (s2/n2) are identical
%  but calculate it *exactly* instead of using randomized permutations.
% 
% i.e. 
%  in group 1 there were s1 successes out of n1 tries,
%  in group 2 there were s2 successes out of n2 tries,
%
%  the returned pvalue is the probability that, if we assign 
%   the (s1 + s2) successes randomly to the two groups, 
%   that the number of successes assigned to group 1 will
%   be as or more extreme as the observed value s1.
%  (where by 'extreme', we (loosely) mean how far it deviates from the
%   expected number of successes in group 1 under the null hypothesis,
%   which would be n1 * ((s1+s2)/(n1+n2)).
%
%  we do this by explicitly calculating the probability of each
%  possible assignment of successes between the two groups.
%
%  This is actually a quite efficient algorithm, if I do say so myself :)
%  The running time and space should both be O(n1+n2).
%  (could super-optimize to reduce space to O(1)! But would have to 
%   call 'log' (n1+n2) times, so it
%   would probably be rather slower)

if nargin < 5
    runmode = 'iterative';
end

switch runmode
    case 'iterative'
        if (s1 == 0 && s2 == 0) || (s1 == n1 && s2 == n2)
            pval = 1.0;
            return;
        end
        
        if (s1 + s2) < (n1 + n2) - (s1 + s2)
            s1 = n1 - s1;
            s2 = n2 - s2;
        end

        % pooled data
        s3 = (s1 + s2);
        n3 = (n1 + n2);

        % minimum/maximum number of successes
        %  that can possibly be assigned to group 1
        min_s1 = max(0,s3-n2);
        max_s1 = min(s3,n1);

        s1hat = min_s1:max_s1;
        s2hat = s3 - s1hat;
        prob_s1hat = nan*ones(size(s1hat));
        log_prob_s1hat = nan*ones(size(s1hat));

        
        % calculate probability that shuffling the successes/failures
        %  would lead to s1hat successes ending up in group 1, and the
        %  remaining (s3 - s1hat) = s2hat successes ending up in group 2:
        %
        %    = nchoosek(n1,s1hat) * nchoosek(n2,s2hat) / nchoosek(n3,s3)
        %
        %    =   n1! / ( s1hat! * (n1 - s1hat)! )
        %      * n2! / ( s2hat! * (n2 - s2hat)! )
        %      * ( s3! * (n3 - s3)! ) / n3!
        %
        % 1. to speed up calculations tremendously, we only calculate 
        %     this equation one time directly, for the smallest
        %     s1hat. After that, we calculate the probability for the
        %     next-biggest s1hat iteratively, using the update:
        %
        %   p(s1hat+1) = p(s1hat) 
        %                * ( (n1 - s1hat) / (s1hat + 1) )
        %                * ( (s2hat) / (n2 - (s2hat-1)) )
        %   
        % 2. to allow for big/small numbers, we do calculations in logspace
        
        
        % -- initial calculation of p(s1hat) for smallest possible s1hat
        
        log_n3_terms_top = log( 1:s3 );
        log_n3_terms_bot = log( (n3-s3+1):n3 );
        log_n3_terms = log_n3_terms_top - log_n3_terms_bot;
        
        i = 1;
        s1h = s1hat(i);
        s2h = s3 - s1h;

        if s1h == 0 || n1-s1h == 0
            log_n1_terms = [];
        else
            log_n1_terms_top = log( (n1-s1h+1):n1 );
            log_n1_terms_bot = log( 1:s1h );
            log_n1_terms = log_n1_terms_top - log_n1_terms_bot;
        end
        if s2h == 0 || n2-s2h == 0
            log_n2_terms = [];
        else
            log_n2_terms_top = log( (n2-s2h+1):n2 );
            log_n2_terms_bot = log( 1:s2h );
            log_n2_terms = log_n2_terms_top - log_n2_terms_bot;
        end
        log_prob_s1hat(i) = sum([log_n1_terms log_n2_terms log_n3_terms]);
        
        
        % -- iterative calculation of the remaining p(s1hat)
        for i = 2:numel(log_prob_s1hat)
            log_prob_s1hat(i) = log_prob_s1hat(i-1) ...
                + log(n1 - s1hat(i-1)) - log(s1hat(i-1)+1) ...
                + log(s2hat(i-1)) - log(n2 - (s2hat(i-1)-1));
        end

        prob_s1hat = exp(log_prob_s1hat);
        
        % -- calculate p-value
        ltprob = sum(prob_s1hat(s1hat < s1));
        eqprob = prob_s1hat(s1hat == s1);
        gtprob = sum(prob_s1hat(s1hat > s1));

        pval = min(1,2*min(ltprob + eqprob,gtprob + eqprob));
        
    case 'iterative_old'
        % Again, correct algorithm but can be very slow
        %  because has to calculate a lot of things redundantly!
        if (s1 + s2) < (n1 + n2) - (s1 + s2)
            s1 = n1 - s1;
            s2 = n2 - s2;
        end

        % pooled data
        s3 = (s1 + s2);
        n3 = (n1 + n2);

        min_s1 = max(0,s3-n2);
        max_s1 = min(s3,n1);

        s1hat = min_s1:max_s1;
        prob_s1hat = nan*ones(size(s1hat));

        log_n3_terms_top = log( 1:s3 );
        log_n3_terms_bot = log( (n3-s3+1):n3 );
        log_n3_terms = log_n3_terms_top - log_n3_terms_bot;
        
        for i = 1:numel(s1hat)
            s1h = s1hat(i);
            s2h = s3 - s1h;

            if (s1h == 0 && s2h == 0) || (n1-s1h == 0 && n2-s2h == 0)
                % not shuffleable
                prob_s1hat(i) = 1;
                continue;
            end

            if s1h == 0 || n1-s1h == 0
                log_n1_terms = [];
            else
                log_n1_terms_top = log( (n1-s1h+1):n1 );
                log_n1_terms_bot = log( 1:s1h );
                log_n1_terms = log_n1_terms_top - log_n1_terms_bot;
            end
            if s2h == 0 || n2-s2h == 0
                log_n2_terms = [];
            else
                log_n2_terms_top = log( (n2-s2h+1):n2 );
                log_n2_terms_bot = log( 1:s2h );
                log_n2_terms = log_n2_terms_top - log_n2_terms_bot;
            end


            % n3_terms is in ascending order.
            %  to 'cancel out' its big entries, we want to
            %  multiply by terms in descending order!
            log_n12_terms = [log_n1_terms log_n2_terms log(ones(1,numel(log_n3_terms)-numel(log_n1_terms)-numel(log_n2_terms)))];
            %log_n12_terms = sort(log_n12_terms,2,'descend');


            log_n_terms = log_n12_terms + log_n3_terms;

            prob_s1hat(i) = exp(sum(log_n_terms));
        end

        ltprob = sum(prob_s1hat(s1hat < s1));
        eqprob = prob_s1hat(s1hat == s1);
        gtprob = sum(prob_s1hat(s1hat > s1));

        pval = min(1,2*min(ltprob + eqprob,gtprob + eqprob));

    case 'basic'
        % CORRECT ALGORITHM, BUT NEEDS TO BE MADE ITERATIVE!
        %  ENDS UP WITH 'INF' VALUES FOR nchoosek AND factorial CALLS!
        % pooled data
        s3 = (s1 + s2);
        n3 = (n1 + n2);

        min_s1 = max(0,s3-n2);
        max_s1 = min(s3,n1);

        s1hat = min_s1:max_s1;
        prob_s1hat = nan*ones(size(s1hat));

        for i = 1:numel(s1hat)
            s1h = s1hat(i);

            prob_s1hat(i) = ( factorial(n3 - s3) * nchoosek(s3,s1h) ) ...
                ./ ( nchoosek(n3,n1) * factorial(n1-s1h) * factorial(n2-(s3-s1h)) );
        end

        ltprob = sum(prob_s1hat(s1hat < s1));
        eqprob = prob_s1hat(s1hat == s1);
        gtprob = sum(prob_s1hat(s1hat > s1));

        pval = min(1,2*min(ltprob + eqprob,gtprob + eqprob));
end
