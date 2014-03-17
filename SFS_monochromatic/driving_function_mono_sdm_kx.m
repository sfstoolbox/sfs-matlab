function D = driving_function_mono_sdm_kx(kx,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_SDM_KX returns the driving signal D for SDM in the kx
%domain
%
%   Usage: D = driving_function_mono_sdm_kx(kx,xs,src,f,[conf])
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
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_SDM_KX(kx,xs,f,src,conf) returns the driving signal for
%   the given secondary sources, desired source type (src), and frequency.
%   To derive the driving signals the spectral division method (SDM) in the kx
%   domain is used.
%
%   see also: plot_sound_field, sound_field_mono_sdm_kx

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargvector(kx);
isargxs(xs);
isargpositivescalar(f);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


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
