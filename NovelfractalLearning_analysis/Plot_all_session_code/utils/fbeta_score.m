function score = fbeta_score(y_predict, y_ture, beta)
% Kaining implemented fbeta score in matlab
% y_predict and y_true binary array with the same size.
assert(all(size(y_predict)==size(y_ture)), "Dimensions of the two array do not match");
if nargin<3 
    beta=2;%default
end

true_pos = sum(y_predict & y_ture);
true_neg = sum(~y_predict & ~y_ture);
false_pos = sum(y_predict & ~y_ture);
false_neg = sum(~y_predict & y_ture);

precision = true_pos/ (true_pos + false_pos);
recall = true_pos / (true_pos + false_neg);

score = (1 + beta^2) * (precision * recall) / (beta^2 * precision + recall);
end