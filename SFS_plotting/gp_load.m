function [x,y,header] = gp_load(file)
%GP_LOAD load x,y from a text file
%   Usage: [x,y,header] = gp_load(file)
%
%   Input parameter:
%       file    - filename of the data file
%
%   Output paramter:
%       x       - x axis values
%       y       - y axis values [vector or matrix]
%       header  - header comment
%
%   GP_LOAD(file) load the values of x and y from a text, that was saved
%   in a Gnuplot compatible format
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate: 2012-02-14 12:06:58 +0100 (Di, 14 Feb 2012) $
% $LastChangedRevision: 634 $
% $LastChangedBy: wierstorf.hagen $


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargfile(file);


%% ===== Computation =====================================================

% FIXME: use fopen etc to read. Handle header line!
% Read the data
if isoctave
    to_be_implemented;
else
    content = textread(file,'','delimiter','\t','commentstyle','shell');
end

x = content(:,1);
y = content(:,2:end);
header = '';
