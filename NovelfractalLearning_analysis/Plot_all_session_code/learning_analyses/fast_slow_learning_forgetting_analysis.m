% This script is to test fast learning - fast forgeting theory
%shuffling_num = 10000;
logical_multiday = cellfun(@(x) ~isempty(x.('FR7410')) | ~isempty(x.('FR7411')), {Neuronlist_all(:).learning})';
logical_multiday = logical_multiday & cellfun(@(x) (numel(x.('learningdate'))==5 || numel(x.('learningdate'))==1 && x.('learningdate')>1), {Neuronlist_all(:).learning})';
logical_multiday = logical_multiday;

novelty_resp_logic = [Neuronlist_all(:).learningforgetinganalysis];

withindaylearningroc = [Neuronlist_all(novelty_resp_logic).withindaylearning_newindex];
acrossdayforgetroc = [Neuronlist_all(novelty_resp_logic).acrossdayforget_newindex];
withindaylearningroc_p = ones(size(withindaylearningroc));% Did not include the p value in this plot.
acrossdayforgetroc_p= ones(size(acrossdayforgetroc));% Did not include the p value in this plot.

neuron_scatterplot([withindaylearningroc' withindaylearningroc_p'], [acrossdayforgetroc', acrossdayforgetroc_p'], ...
    'Within day learning index', 'Across day forgetting index', [-1.1,1.1], [-1.1,1.1], plotpath, 'Correlation of learning and forgetting.pdf');
