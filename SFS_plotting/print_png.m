function print_png(filename)
%PRINT_PNG creates a png file from the last plot command
%   Usage: print_png(filename)
%
%   Input options:
%       filename    - name of the output png file
%
%   PRINT_PNG(filename) creates a png version from the last plot command in the
%   in the given file. The size of the png will be 500x375 pixel. The line width
%   etc. of the axis is arranged in order to create a nice looking plot.
%
%   see also: print_eps

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
if ~ischar(filename)
    error('%s: filename has to be a string!',upper(mfilename));
end


%% ===== Plotting ========================================================
% Save old values in order to restore them
linewidth = get(gca,'LineWidth');
ticklength = get(gca,'TickLength');
position = get(gca,'Position');
set(gca,'LineWidth',2);
set(gca,'TickLength',[0.005 0.005]);
% FIXME: check if the following also doesn't work under Matlab!
% Update Position by hand to aviod clipping of labels
%sc = 1.2;
%set(gca,'Position',...
%    [position(1)*sc,position(2)*sc,position(3)/sc,position(4)/sc]);
print(filename,'-dpng','-r75');
set(gca,'LineWidth',linewidth);
set(gca,'TickLength',ticklength);
set(gca,'Position',position);
