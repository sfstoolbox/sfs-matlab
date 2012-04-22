function h=figsize(x,y);
% function h=figsize(x,y);
% Creates a new figure with x and y dimensions in cm

% FIXME: check if function can be removed

if nargin==0
  x=10;
  y=10;
end
if nargin==1
  y=x;
end
% h=figure;
% zoom on;
% colormap(1-gray);
set(gcf,'Position',[20,45,44*x,44*y]);
x=x+2;
y=y+2;
%set(gca,'FontName','Times');
%set(gca,'FontSize',8);
% set(gca,'Position',[1/x,1/y,1-2/x,1-2/y]);
set(gcf,'PaperUnits','centimeters');
% FIXME: the following doesn't work with Octave (is this needed at all?)
%set(gcf,'Papertype','A4')
tmp=get(gcf,'Papersize');
set(gcf,'PaperPosition',[(tmp(1)-x)/2,(tmp(2)-y)/2,x,y]);
set(gcf,'Color',[1,1,1]);

