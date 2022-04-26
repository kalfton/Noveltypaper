function [Hist_raster, ax] = Raster_PSTH_plot(sptimes, time_start, Hist_template, ax)
% sptimes: cell array, each cell contains the spike times (unit: second) of a
% trial. 
% time_start: an array of start times (unit s) to allign each trial on. 
% Hist_template: an array of times (unit ms) relative to the start time,
% example input: -200:2500, mean plot the Raster from -200ms to 2500ms, relative to the start times.
% Kaining modified from Ilya's code
addpath('Y:\Kaining\HELPER_GENERAL\');

if numel(sptimes)~=numel(time_start)
    error('The input size is not the same.')
end

if ~exist('Hist_template', 'var') || isempty(Hist_template)
    Hist_template = -200: 2500;
end
Hist_raster = [];
Hist_sdf = [];
for x=1:numel(time_start)
    Hist_raw = histcounts((sptimes{x}(1:end)-time_start(x))*1000, Hist_template);
    Hist_raster = [Hist_raster; Hist_raw];
    
    Smoothing=2; %1 - box car; 2- gaus; 3 - epsp (causal); specs of kernel are defined in Smooth_Histogram function
    sdf=Smooth_Histogram(Hist_raw,Smoothing);
    Hist_sdf = [Hist_sdf;sdf];
    %clear temp spk x
end

%% plot
if ~exist('ax', 'var') || isempty(ax)
    figure;
    ax = subplot(1,1,1);
end

h1=area(nanmean(Hist_sdf,1));
h1.FaceColor = 'red';

xt = sum(Hist_raster>0,2); % number of spikes in a trial

rastList = nan(size(Hist_raster,1),max(xt));
for tq=1:size(Hist_raster,1)
    temptq=find(Hist_raster(tq,:)>0);
    rastList(tq,1:length(temptq))=temptq;
end

rasIntv=1;
LWidth=1;
LColor='k';
maxY_rast=20;

for line = 1:size(rastList,1)
    hold on
    curY_rast = maxY_rast+rasIntv*line;
    nonnanind = ~isnan(rastList(line,:));
    plot([rastList(line,nonnanind); rastList(line,nonnanind)],...
        [(curY_rast+rasIntv/2)*ones(1,sum(nonnanind)); ...
        (curY_rast-rasIntv/2)*ones(1,sum(nonnanind))],'-','linewidth',LWidth,'color',LColor )
end

end
