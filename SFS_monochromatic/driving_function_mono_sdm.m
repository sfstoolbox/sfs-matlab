function D = driving_function_mono_sdm(x0,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_SDM returns the driving signal D for SDM
%
%   Usage: D = driving_function_mono_sdm(x0,xs,src,f,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx7]
%       xs          - position of virtual source or direction of plane
%                     wave / m [1x3]
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_SDM(x0,xs,f,src,conf) returns the driving signal for
%   the given secondary sources, desired source type (src), and frequency.
%   To derive the driving signals the spectral division method (SDM) is used.
%
%   See also: plot_sound_field, sound_field_mono_sdm, driving_function_imp_sdm

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
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargxs(xs);
isargpositivescalar(f);
isargchar(src);
isargstruct(conf);


%% ===== Computation ====================================================

% Calculate the driving function in time-frequency domain

% Secondary source positions and directions
nx0 = x0(:,4:6);
x0 = x0(:,1:3);

% Source position
xs = repmat(xs(1:3),[size(x0,1) 1]);

% Get driving signals
if strcmp('pw',src)
    % === Plane wave =====================================================
    % Direction of plane wave
    nk = bsxfun(@rdivide,xs,vector_norm(xs,2));
    % Driving signal
    D = driving_function_mono_sdm_pw(x0,nk,f,conf);

elseif strcmp('ps',src)
    % === Point source ===================================================
    % Driving Signal
    D = driving_function_mono_sdm_ps(x0,xs,f,conf);

elseif strcmp('ls',src)
    % === Line source ====================================================
    % Driving signal
    D = driving_function_mono_sdm_ls(x0,xs,f,conf);

elseif strcmp('fs',src)
    % === Focused source =================================================
    % Driving Signal
    D = driving_function_mono_sdm_fs(x0,xs,f,conf);

else
    error('%s: %s is not a known source type.',upper(mfilename),src);
end
