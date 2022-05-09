function [h] = plot_raster(x,y,plotstyle,varargin)
% [h] = plot_raster(x,y,plotstyle, ...)
%
% plots a spike raster using various styles:
%
% plot_raster(x,y,'point',mark,...)
% => plot(x,y,'linestyle','none','marker',mark,...)
%  (default marker: '.')
%
% plot_raster(x,y,'segment',spike_height,...)
% => plot(x,y,...), where
%  x = [x(1) x(1) nan x(2) x(2) nan ...]
%  y = [y(1) y(1) nan y(2) y(2) nan ...] + spike_height*[+.5 -.5 nan +.5 -.5 nan ...]
% (default spike_height: 1)
%
% plot_raster(x,y,'patch',spike_width,spike_height,spike_color,...)
% => 
% (default spike_width 1, spike_height 1, spike_color 'k')
%
% plot_raster(x,y,'image',x_smooth_style,y_smooth_style,spike_width,spike_color)
% => creates image version and plots the image
%
%
% outputs:
%  h - handle to the lineseries

assert(numel(x) == numel(y));

switch plotstyle
    case 'point'
        if numel(varargin) >= 1
            mark = varargin{1};
            varargin = varargin(2:end);
        else
            mark = '.';
        end;
        h=plot(x,y,'linestyle','none','marker',mark,varargin{:});
        
    case 'segment'
        if numel(varargin) >= 1
            spike_height = varargin{1};
            varargin = varargin(2:end);
        else
            spike_height = 1;
        end;
        
        xnew = nans(1,numel(x)*3);
        ynew = nans(1,numel(y)*3);
        
        cur_i = 1 + ((1:numel(x)) - 1)*3;
        xnew(cur_i) = x;
        ynew(cur_i) = y + (spike_height/2);
        
        cur_i = 2 + ((1:numel(x)) - 1)*3;
        xnew(cur_i) = x;
        ynew(cur_i) = y + (spike_height/2);
        
        for i = 1:numel(x)
            xnew(1 + (i-1)*3) = x(i);
            xnew(2 + (i-1)*3) = x(i);
            ynew(1 + (i-1)*3) = y(i)+(spike_height/2);
            ynew(2 + (i-1)*3) = y(i)-(spike_height/2);
        end;
        
        h=plot(xnew,ynew,'-',varargin{:});
        
    case 'patch'
        if numel(varargin) >= 1
            spike_width = varargin{1};
            varargin = varargin(2:end);
        else
            spike_width = 1;
        end;
        if numel(varargin) >= 1
            spike_height = varargin{1};
            varargin = varargin(2:end);
        else
            spike_height = 1;
        end;
        if numel(varargin) >= 1
            spike_color = varargin{1};
            varargin = varargin(2:end);
        else
            spike_color = 'k';
        end;
        
        h = nans(size(x));
        for i = 1:numel(x)
            h(i) = plot_rectpatch(...
                x(i) + spike_width*0.5*[-1 1], ...
                y(i) + spike_height*0.5*[-1 1], ...
                spike_color, varargin{:});
        end;
        
    case 'image'
        % plot_raster(x,y,'image',x_smooth_style,y_smooth_style,spike_width,spike_color)
        if numel(varargin) >= 1
            x_smooth = varargin{1};
            varargin = varargin(2:end);
        else
            x_smooth = 'none';
        end;
        if numel(varargin) >= 1
            y_smooth = varargin{1};
            varargin = varargin(2:end);
        else
            y_smooth = 'none';
        end;
        if numel(varargin) >= 1
            spike_width = varargin{1};
            varargin = varargin(2:end);
        else
            spike_width = 1;
        end;
        if numel(varargin) >= 1
            spike_color = varargin{1};
            varargin = varargin(2:end);
        else
            spike_color = 'k';
        end;
        
        x = x(:);
        y = y(:);
        
        xrange = [(min(x)-spike_width) (max(x)+spike_width)];
        yrange = [(min(y)-1) (max(y)+1)];
        xbin = xrange(1):spike_width:xrange(2);
        ybin = yrange(1):yrange(2);
        hm = hist2(y,x,ybin,xbin);
        
        for ss = 1:2
            switch ss
                case 1
                    % smooth along X dimenson
                    smooth = x_smooth;
                case 2
                    % smooth along Y dimenson
                    smooth = y_smooth;
                    hm = hm';
            end;
            
            if ~isequal(smooth,'none')
                smooth_style = smooth{1};
                smooth_param = smooth{2};
                switch smooth_style
                    case 'runningaverage'
                        if iscell(smooth_param)
                            hm = runningaverage(hm,smooth_param{1});
                            hm = hm .* smooth_param{2};
                        else
                            hm = runningaverage(hm,smooth_param);
                            hm = hm .* smooth_param;
                        end;
                    case 'gauss'
                        if iscell(smooth_param)
                            [hm,used_filter] = efilter(hm,'gauss',smooth_param{1});
                            hm = hm .* smooth_param{2};
                        else
                            [hm,used_filter] = efilter(hm,'gauss',smooth_param);
                            hm = hm ./ max(used_filter);
                        end;
                    case 'rebin'
                        if iscell(smooth_param)
                            len_to_rebin = floor(size(hm,2)./smooth_param{1}).*smooth_param{1};
                            hm = [ ...
                                rebin(hm(:,1:len_to_rebin),smooth_param{1},@mean,false,false) ...
                                zeros(size(hm,1),size(hm,2)-len_to_rebin)];
                            hm = hm .* smooth_param{2};
                        else
                            len_to_rebin = floor(size(hm,2)./smooth_param).*smooth_param;
                            hm = [ ...
                                rebin(hm(:,1:len_to_rebin),smooth_param,@max,false,false) ...
                                zeros(size(hm,1),size(hm,2)-len_to_rebin)];
                        end;
                    otherwise
                        error('unknown x smoothing style');
                end;
            end;
        
            switch ss
                case 2
                    hm = hm';
            end;
        end;
        
        
        h=image(xrange,yrange,colormapify(hm,[0 1],'w',interpcolor('w',spike_color,0.5),spike_color));
        
    otherwise
        error('unknown plotstyle');
end;
