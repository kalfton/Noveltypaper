% This script is to test fast learning - fast forgeting theory
%shuffling_num = 10000;
logical_multiday = cellfun(@(x) ~isempty(x.('FR7410')) | ~isempty(x.('FR7411')), {Neuronlist_good(:).learning})';
logical_multiday = logical_multiday & cellfun(@(x) (numel(x.('learningdate'))==5 || numel(x.('learningdate'))==1 && x.('learningdate')>1), {Neuronlist_good(:).learning})';
logical_multiday = logical_multiday;

novelty_resp_logic = [Neuronlist_good(:).learningforgetinganalysis];

withindaylearningroc = [Neuronlist_good(novelty_resp_logic).withindaylearning_newindex];
acrossdayforgetroc = [Neuronlist_good(novelty_resp_logic).acrossdayforget_newindex];%.*withindaylearningroc;
withindaylearningroc_p = ones(size(withindaylearningroc));%[Neuronlist_good(novelty_resp_logic).withindaylearningroc_p];
acrossdayforgetroc_p= ones(size(acrossdayforgetroc));%[Neuronlist_good(novelty_resp_logic).acrossdayforgetroc_p];

scatterplot_goodlooking([withindaylearningroc' withindaylearningroc_p'], [acrossdayforgetroc', acrossdayforgetroc_p'], ...
    'Within day learning index', 'Across day forgetting index', [-1.1,1.1], [-1.1,1.1], plotpath);
