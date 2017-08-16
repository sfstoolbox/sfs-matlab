function [d,x0,xv,idx,delay_offset] = driving_function_imp_localwfs_vss(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_LOCALWFS_VSS returns the driving signal d for local WFS
%using focused sources as virtual secondary sources
%
%   Usage: [d,x0,xv,idx,delay_offset] = driving_function_imp_localwfs_vss(x0,xs,src,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx7]
%       xs          - position of virtual source or direction of plane
%                     wave / m [1x3]
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'ls' - line source
%                         'fs' - focused source
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       d             - driving signals [mxn]
%       x0            - position, direction, and weights of the real secondary
%                       sources / m [nx7]
%       xv            - position, direction, and weights of the virtual 
%                       secondary sources / m [mx7]
%       idx           - index of the selected sources from the original x0
%                       matrix [mx1]
%       delay_offset  - additional added delay, so you can correct it
%
%   References:
%       S. Spors (2010) - "Local Sound Field Synthesis by Virtual Secondary
%                          Sources", 40th AES
%
%   See also: plot_sound_field, sound_field_mono_wfs

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargxs(xs);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ==================================================
virtualconf = conf;
virtualconf.usetapwin = conf.localwfs_vss.usetapwin;
virtualconf.tapwinlen = conf.localwfs_vss.tapwinlen;
virtualconf.secondary_sources.size = conf.localwfs_vss.size;
virtualconf.secondary_sources.center = conf.localwfs_vss.center;
virtualconf.secondary_sources.geometry = conf.localwfs_vss.geometry;
virtualconf.secondary_sources.number = conf.localwfs_vss.number;
virtualconf.wfs = conf.localwfs_vss.wfs;
virtualconf.nfchoa = conf.localwfs_vss.nfchoa;
virtualconf.driving_functions = conf.localwfs_vss.driving_functions;
method = conf.localwfs_vss.method;

%% ===== Computation ====================================================
if strcmp('fs',src)
    error(['%s: %s is not a supported method source type! Try to use a', ...
        ' point source, if the source is inside the secondary source array', ...
        ' but not inside the virtual secondary source array'], ...
        upper(mfilename),src);
end

% Determine driving functions of virtual array with different sfs methods
switch method
case 'wfs'
    % === Wave Field Synthesis ===
    % Create virtual source array
    xv = virtual_secondary_source_positions(x0,xs,src,conf);
    % Secondary_source_selection
    xv = secondary_source_selection(xv, xs, src);
    % Optional tapering
    xv = secondary_source_tapering(xv,virtualconf);
    % Driving functions for virtual source array
    [dv,~,~,delay_offset] = driving_function_imp_wfs(xv,xs,src,virtualconf);   
case 'nfchoa'
    % === Near-Field-Compensated Higher Order Ambisonics ===
    % Create virtual source array
    xv = secondary_source_positions(virtualconf);
    % Driving functions for virtual source array 
    [dv,~,delay_offset] = driving_function_imp_nfchoa(xv,xs,src,virtualconf);
otherwise
    error('%s: %s is not a supported method for localsfs!',upper(mfilename),...
      method);
end

% Select secondary sources
[x0,idx] = secondary_source_selection(x0,xv(:,1:6),'vss');
% Driving functions for real source array
[d,delay_vss] = driving_function_imp_wfs_vss(x0,xv,'fs',dv,conf);
% add delay
delay_offset = delay_offset + delay_vss;
