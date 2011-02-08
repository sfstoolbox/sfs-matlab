function [X,Y] = setting_xy_ranges(X,Y,conf)
%SETTING_XY_RANGES sets the ranges of the X and Y axes in the wave field
%simulations
%   Usage: [X,Y] = setting_xy_ranges(X,Y,conf)
%          [X,Y] = setting_xy_ranges(X,Y)
%
%   Input options:
%       X,Y     -   X,Y values. Could be in the form of [xmin,xmax],[ymin,ymax] 
%                   or single values X,Y.
%       conf    -   optional configuration struct (see SFS_config)
%
%   Output options:
%       X,Y -   x- and y-range in the form of [xmin,xmax],[ymin,ymax]
%
%   SETTING_XY_RANGES(X,Y) returns two verctors containing the x- and y-range
%   needed for the wave field simulations. If X and Y are only single values,
%   then the range is calculated depending on the used array (see conf.array).
%   xmin will be always -X/2 and xmax X/2, but ymin depends on the array form
%   and will be -0.1 for a linear array and -Y/2 for a circular array.
%
%   see also: wave_field_time_domain, wave_field_monochromatic_wfs_25d

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
isargvector(X,Y);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================

array = conf.array;
lssize = conf.plot.lssize;


%% ===== Computation =====================================================

% Setting the x- and y-axis
if isscalar(X)
    tmp = X;
    X(1) = -tmp/2;
    X(2) = tmp/2;
elseif length(X)>2
    error('%s: X has to be a single value or [xmin xmax]!',upper(mfilename));
end
if isscalar(Y) && strcmp('linear',array)
    tmp = Y;
    Y(1) = -lssize;
    Y(2) = tmp;
elseif isscalar(Y)
    tmp = Y;
    Y(1) = -tmp/2;
    Y(2) = tmp/2;
elseif length(Y)>2
    error('%s: Y has to be a single value or [ymin ymax]!',upper(mfilename));
end
