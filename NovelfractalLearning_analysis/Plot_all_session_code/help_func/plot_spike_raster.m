%******************* make a rastergram of the actual spikes on each trial
function plot_spike_raster(spmat,xstart,xend,TP)
%*********** spmat is a 1xC cell with spike matrices of conditions 
%***********     for which you wish to smooth the rate and plot
%***********     the mean firing rate, each matrix is (trials x time)
%*********** xstart - start of interval to plot
%*********** xend - end of interval to plot
%*********** TP - time sampling, defaults to TP = 1


if (isempty(TP))
    TP = 1;
    ixstart = xstart;
    ixend = xend;
else
    ixstart = 1;
    ixend = size(spmat{1,1},2);
end

col = 'rbg';
colo = [[1,1,1];[0.5,0,0];[0,0,0.5];[0,0.5,0]];
colormap(colo);


totspike = [];
for cubo = 1:size(spmat,2)
    totspike = [totspike; (cubo*spmat{1,cubo}(:,ixstart:ixend)) ];
    totspike = [totspike; zeros(2,size(totspike,2))];
end
totspike = totspike + 1;
imagesc(xstart:(xstart+(TP*(ixend-ixstart))),...
          (0:(size(totspike,1))+1),totspike); hold on;
%*************
it = 0.5;

for cubo = 1:size(spmat,2)
    plot([xstart,xend],[it,it],[col(cubo),'-']);
    it = it + 1 + size(spmat{1,cubo},1);
    plot([xstart,xend],[it,it],[col(cubo),'-']);
    it = it + 1;
end

return;

