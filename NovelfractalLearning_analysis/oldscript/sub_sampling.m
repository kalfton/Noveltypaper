function select_ind = sub_sampling(data, gaussian_mean, gaussian_std, varargin)
% input a sequence of indices of neuron, the goal of this function is to
% subsample the neurons, such that the indices of the subsample has the
% target probability distribution.
n_div = 100;
min_thr = 1e-5;
data_range = [min(data), max(data)];
% gaussian_mean = 0.4;
% gaussian_std = 0.03;
bins = data_range(1):gaussian_std/n_div:data_range(2);
bincenter = bins(1:end-1)+gaussian_std/n_div/2;

datahist = histcounts(data, bins);
%smooth the histogram using gaussian kernel
datahist = conv(datahist,gausswin(10*n_div));
datahist = datahist(5*n_div:end-5*n_div);


% taget distribution:
targetpdf = normpdf(bincenter, gaussian_mean, gaussian_std)';


% if the numbers are too small, set them to zero.
targetpdf(datahist<min_thr) = 0;
datahist(datahist<min_thr) = 0;
% find the points that datahist(i) = 0 while targetpdf(i) \neq 0, if exist
% at least a point return error.
if ~isempty(find(datahist==0 & targetpdf~=0, 1))|| gaussian_mean-data_range(1)<2.5*gaussian_std || data_range(2)-gaussian_mean<2.5*gaussian_std
    error('the target distribution is out of the range of the data');
end

% make the two arrays in the same dimension.
if size(targetpdf,1) ~=size(datahist,1)
    targetpdf = targetpdf';
end
% get the accept/reject rate at each point and rescale it
acceptrate = targetpdf./datahist;
acceptrate(datahist<min_thr)=0;
acceptrate(isnan(acceptrate))=0;
%


% normalize the rate
acceptrate = acceptrate/max(acceptrate);

% subsampling the data
select_logical = zeros(size(data));
for i = 1:numel(data)
    [~, ix] = min( abs( data(i)-bincenter ) );
    rate = acceptrate(ix);
    if rand<rate
        select_logical(i) = 1;
    end
end
select_ind = find(select_logical);






end

% data = randn(1000,1);
% figure;
% nsubplot(110,110, 10:50, 10:50);
% histogram(data);
% xlim([-5,5]);
% 
% temp_ind = sub_sampling_fun(data, 2, 0.1);
% sub_data = data(temp_ind);
% 
% nsubplot(110,110, 50+(10:50), 10:50);
% histogram(sub_data);
% xlim([-5,5]);

