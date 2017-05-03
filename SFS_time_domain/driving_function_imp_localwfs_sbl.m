function [d,delay_offset] = driving_function_imp_localwfs_sbl(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_LOCALWFS_SBL returns the driving signal for local WFS
%using spatial bandwidth limitation
%
%   Usage: [d,delay_offset] = driving_function_imp_localwfs_sbl(x0,xs,src,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx7]
%       xs          - position of virtual source or direction of plane
%                     wave / m [1x3]
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       d               - driving signals [mxn]
%       delay_offset    - additional added delay, so you can correct it
%
%   See also: sound_field_imp, sound_field_imp_localwfs_sbl

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


%% ===== Configuration ========================================================
t0 = conf.t0;


%% ===== Computation ==========================================================

% Needed to time-align lf and hf part of driving function for point source
if ~strcmp(t0, 'source')
    error('%s: conf.t0 (%s) other than "source" is not supported', ...
        upper(mfilename),t0);
end

switch src
case 'ps'
    [d,delay_offset] = driving_function_imp_localwfs_sbl_ps(x0,xs,conf);
case 'pw'
    xs = xs./norm(xs);
    [d,delay_offset] = driving_function_imp_localwfs_sbl_pw(x0,xs,conf);
otherwise
    error('%s: %s is not a known source type.',upper(mfilename),src);
end
