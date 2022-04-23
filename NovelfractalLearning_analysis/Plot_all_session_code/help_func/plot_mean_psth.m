%****************** function to make a nice raster plot ****************
function TRIALS = plot_mean_psth(spmat,smooth_window,xstart,xend,TP)
%*********** spmat is a 1xC cell with spike matrices of conditions 
%***********     for which you wish to smooth the rate and plot
%***********     the mean firing rate, each matrix is (trials x time)
%*********** smooth_window - sigma for gaussian window smoothing
%*********** xstart - start of interval to plot
%*********** xend - end of interval to plot
%*********** TP - time precision, default is 1 ms

numconds = size(spmat,2);
if (numconds==2)
    colo = [[1,0,0];[0,0,1]];
else
    colo = [[1,0,0];[0,0,1];[0,1,0];[1,1,0]];
end

if (isempty(TP)) | (TP == 1)
    TP = 1;
    ixstart = xstart;
    ixend = xend;
else
    ixstart = 1;
    ixend = size(spmat{1,1},2);
end

maxo = 0;
for k = 1:numconds
  spud = spmat{1,k};  
  numtrials(k) = size(spud,1);
  smorate = compute_gauss_smooth(sum( spud(1:numtrials(k),ixstart:ixend))/....
                      numtrials(k),(smooth_window/TP))*1000;
%  H = plot(xstart:TP:(xstart+(TP*(ixend-ixstart))),smorate,'k-'); hold on;
%  set(H,'Color',colo(k,:));
  
  if (max(smorate) > maxo)
      maxo = max(smorate);
  end
end

numtrials_(k) = size(spud,1);
TRIALS=[];
for k = 1:numtrials_
  temp=spud(k,:);
  numtrials(k) = size(temp,1);
  smorate = compute_gauss_smooth(( temp(1:numtrials(k),ixstart:ixend))/....
                      numtrials(k),(smooth_window/TP))*1000;
  TRIALS=[TRIALS; smorate']; 
  clear smorate
end


return;


