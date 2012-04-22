% appends a colorbar with title to the current figure
% S.Spors / 24.9.2004

% FIXME: this function is not for the public yet
function []=tcolorbar(titlel,titlet)

hcbar = colorbar;

h=get(hcbar,'YLabel');
set(h,'string',titlel);
pos=get(h,'Position');
set(h,'Position',pos+[0.3 0 0]);

if nargin==2
  set(get(hcbar,'Title'),'String',titlet);
end
