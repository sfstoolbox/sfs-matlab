function fontsize(fsize)
%FONTSIZE sets the size of the figure fonts
%   Usage: fontsize(fsize)
%
%   Input options:
%       fsize   - font size
%
%   FONTSIZE(fsize) sets the fonts of the active figure to the given size.
%
%   See also: GraphDefaults, fontname

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(fsize) || ~isscalar(fsize) || fsize<=0
    error('%s: fsize has to be a positive scalar!',upper(mfilename));
end


%% ===== Apply settings ==================================================
% Get handle for active figure
h=gca;
% Set the font size for the figure
set(h,'FontSize',fsize);
temp=get(h,'Xlabel');
xlabel(get(temp,'String'));
set(temp,'FontSize',fsize);
temp=get(h,'Ylabel');
ylabel(get(temp,'String'));
set(temp,'FontSize',fsize);
temp=get(h,'Title');
title(get(temp,'String'));
set(temp,'FontSize',fsize);
