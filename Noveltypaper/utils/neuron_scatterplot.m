function neuron_scatterplot(x_data, y_data, x_name, y_name, x_limit, y_limit, plotpath, plotname)
% x_data: n*2 array, first column index, second column p value
% y_data: n*2 array, first column index, second column p value
% x_name: the name of x data
% y_name: the name of y data

HISTWIDTH=0.05; 
edges = linspace(x_limit(1), x_limit(2), 21); % Create 20 bins.
THR=0.01; %stat threshold


figure;

nsubplot(5, 5, 1, 3); set(gca,'ticklength',4*get(gca,'ticklength'));

datatoplot=x_data;
totalleng = size(datatoplot, 1);

h1 = histogram(datatoplot(:,1), 'EdgeColor','k','FaceColor','k','BinEdges',edges); hold on
h1.BinWidth = HISTWIDTH;

h2 = histogram(datatoplot(find(datatoplot(:,2)<THR),1), 'EdgeColor','k','FaceColor','w','BinEdges',edges); hold on
h2.BinWidth = HISTWIDTH;

yl =  ylim; YL=round(yl(2));
xlim(x_limit);
ylim([0 YL]);
XL =  xlim; XL=round(XL(2));

scatter(nanmean((datatoplot(:,1))),YL,80,'r','v','filled')
line([nanmean(datatoplot(:,1)) nanmean(datatoplot(:,1))], [0 YL] , 'Color', [0 1 0 ],'LineWidth',1) ;

%%
x=1; pr=((YL)./5)*x; text(XL,pr,['n= ' mat2str(length(datatoplot)) '/' mat2str(totalleng)])

[p,h]=signrank(datatoplot(:,1),0);
x=2; pr=((YL)./5)*x; text(XL,pr,['p = ' mat2str(p)]);
line([0 0], [0 YL] , 'Color', [1 0 0 ],'LineWidth',2) ;

ylabel(['Count'])
axis square;
 

    
nsubplot(5, 5, 3, 1); set(gca,'ticklength',4*get(gca,'ticklength'));

datatoplot=y_data;
totalleng = size(datatoplot, 1);

h1 = histogram(datatoplot(:,1), 'EdgeColor','k','FaceColor','k','BinEdges',edges); hold on
h1.BinWidth = HISTWIDTH;

h2 = histogram(datatoplot(find(datatoplot(:,2)<THR),1), 'EdgeColor','k','FaceColor','w','BinEdges',edges); hold on
h2.BinWidth = HISTWIDTH;

yl =  ylim; YL=round(yl(2));
xlim(y_limit);
ylim([0 YL]);
XL =  xlim; XL=round(XL(2));

line([nanmean(datatoplot(:,1)) nanmean(datatoplot(:,1))], [0 YL] , 'Color', [0 1 0 ],'LineWidth',1) ;
line([0 0], [0 YL] , 'Color', [1 0 0 ],'LineWidth',2) ;

ylabel(['Count'])

set(gca,'view',[-90 90])
axis square;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

nsubplot(5, 5, 4, 1); set(gca,'ticklength',4*get(gca,'ticklength'));

x=1; pr=((YL)./5)*x; text(XL,pr,['n= ' mat2str(length(datatoplot)) '/' mat2str(totalleng)])

[p,h]=signrank(datatoplot(:,1),0);
x=2; pr=((YL)./5)*x; text(XL,pr,['p = ' mat2str(p)])

xlim(y_limit);
ylim([0 YL]);

axis off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
nsubplot(5, 5, 3, 3); set(gca,'ticklength',4*get(gca,'ticklength'));
axis square;
%ScatterIndicesWithColors(x_data,y_data);
show_cb =1;
ScatterIndicesheatmap(x_data,y_data, x_limit, show_cb);
xlabel(x_name);
ylabel(y_name);
xlim(x_limit);
ylim(y_limit);


set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figuren and orientation
if exist('plotpath', 'var')
    print(gcf,'-dpdf', '-painters',fullfile(plotpath, plotname));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function ScatterIndicesheatmap(datax,datay,lim, show_cb)
    %lim = [-1, 1];
    nbin = 25;
    binedge = linspace(lim(1),lim(2),nbin+1);
    binmid = 0.5*binedge(1:(end-1)) + 0.5*binedge(2:end);
    
    [n,xedge,yedge] = histcounts2(datax(:,1),datay(:,1),binedge,binedge);
    
    n = n'; % transpose so that x=columns and y=rows
    n_origin = n;
    
    n = n ./ nansum(n(:));
    n = log(n);
    minimum_n = min(n(n>-inf));
    n(n==-inf)=minimum_n-0.5;
    colormap(flipud(gray));
    
    imagesc(binmid,binmid,n);
    n_contour_levels = 6;
    contour(binmid,binmid,n,n_contour_levels,'linewidth',2);
    
    if show_cb
        cb = colorbar;
        cb.Ticks = log([1,10,100,1000]./nansum(n_origin(:)));
        cb.TickLabels = round(exp(cb.Ticks)*nansum(n_origin(:)));
    end
    
    try
    
    line(lim,[0 0],'color',[0.3 0.3 0.3],'LineWidth',1);
    line([0 0],lim,'color',[0.3 0.3 0.3],'LineWidth',1);
    axis([-0 1 0 1])
    
    xticks([-1 :0.5: 1])
    yticks([-1 :0.5: 1])
    axis square
    
    
    [rho,pval] = corr(datax(:,1),datay(:,1), 'Type', 'Spearman');
    
    maxval=ylim;
    text(1,0.7,['p = ' mat2str(pval,3)],'FontSize',11)
    text(1,0.6,['rho = ' mat2str(rho,3)],'FontSize',11)
    text(maxval(1), maxval(2), [' n = ' mat2str(length(datax)) ])
    
    index=(datax(:,2)<0.05 | datay(:,2)<0.05);
    [rho,pval] = corr(datax(index,1),datay(index,1), 'Type', 'Spearman');
    
    text(1,1,['p_sig = ' mat2str(pval,3)],'FontSize',11)
    text(1,0.8,['rho_sig = ' mat2str(rho,3)],'FontSize',11)
end
