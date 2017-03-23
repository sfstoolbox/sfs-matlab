function D = driving_function_mono_sdm_kx(kx,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_SDM_KX returns the driving signal D for SDM in the kx
%domain
%
%   Usage: D = driving_function_mono_sdm_kx(kx,xs,src,f,conf)
%
%   Input parameters:
%       kx          - kx dimension [nx1]
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
%   DRIVING_FUNCTION_MONO_SDM_KX(kx,xs,f,src,conf) returns the driving signal for
%   the given secondary sources, desired source type (src), and frequency.
%   To derive the driving signals the spectral division method (SDM) in the kx
%   domain is used.
%
%   See also: plot_sound_field, sound_field_mono_sdm_kx

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
isargvector(kx);
isargxs(xs);
isargpositivescalar(f);
isargchar(src);
isargstruct(conf);


%% ===== Computation ====================================================

% Calculate the driving function in time-frequency domain

% Get driving signals
if strcmp('pw',src)
    % === Plane wave =====================================================
    % Direction of plane wave
    nk = xs / norm(xs);
    % Driving signal
    D = driving_function_mono_sdm_kx_pw(kx,nk,f,conf);

elseif strcmp('ps',src)
    % === Point source ===================================================
    % Driving Signal
    D = driving_function_mono_sdm_kx_ps(kx,xs,f,conf);

elseif strcmp('ls',src)
    % === Line source ====================================================
    % Driving signal
    D = driving_function_mono_sdm_kx_ls(kx,xs,f,conf);

elseif strcmp('fs',src)
    % === Focused source =================================================
    % Driving Signal
    D = driving_function_mono_sdm_kx_fs(kx,xs,f,conf);

else
    error('%s: %s is not a known source type.',upper(mfilename),src);
end
