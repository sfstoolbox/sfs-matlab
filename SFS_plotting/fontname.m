function fontname(fname)
%FONTSIZE sets the figures font
%   Usage: fontname(fname)
%
%   Input options:
%       fname   - font name
%
%   FONTNAME(fname) sets the font of the active figure to the given type.
%
%   See also: GraphDefaults, fontsize

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));

if ~ischar(fname)
    error('%s: fname has to be a string!',upper(mfilename));
end


%% ===== Apply settings ==================================================
% Get handle for active figure
h=gca;
% Set the font type 
set(h,'FontName',fname);
temp=get(h,'Xlabel');
xlabel(get(temp,'String'));
set(temp,'FontName',fname);
temp=get(h,'Ylabel');
ylabel(get(temp,'String'));
set(temp,'FontName',fname);
temp=get(h,'Title');
title(get(temp,'String'));
set(temp,'FontName',fname);
