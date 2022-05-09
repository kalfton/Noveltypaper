function [h] = figuren(varargin)
% [h] = figuren(varargin)
% same as 'figure', but applies my favorite settings using 'nsubplot'!

h = figure(varargin{:});set(gcf,'color',[1 1 1])
nsubplot;