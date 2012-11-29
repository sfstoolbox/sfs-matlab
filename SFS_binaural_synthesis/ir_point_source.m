function ir_ps = ir_point_source(X,phi,xs,irs,conf)
%IR_POINT_SOURCE Generate a IR for a point source
%
%   Usage: ir_ps = ir_point_source(X,phi,xs,irs,[conf])
%
%   Input parameters:
%       X       - listener position (m)
%       phi     - listener direction [head orientation] (rad)
%       xs      - source position (m)
%       irs     - IR data set for the second sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       ir_ps   - Impulse response (nx2 matrix)
%
%   IR_POINT_SOURCE(X,phi,xs,irs,conf) calculates a impulse response for a
%   reference source (single loudspeaker) at position xs and a listener
%   located at X and looking into direction phi. Whereby at phi = 0 the
%   listener is looking in the direction of the x-axis, like the angle is
%   normally defined in Mathematics.
%
%   Geometry:
%
%                                 y-axis
%                                   ^
%                                   |
%                                   |
%                                   |
%                                   |    listener
%                                   |       O X, phi=-pi/2
%                                   |       |
%               source              |
%                 o xs              |
%                                   |
%                                   |
%       ----------------------------|---------------------------> x-axis
%
%
% see also: SFS_config, ir_wfs_25d, auralize_ir, brs_point_source

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ====================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
[X,xs] = position_vector(X,xs);
isargscalar(phi);
check_irs(irs);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
fs = conf.fs;                  % sampling frequency
c = conf.c;                    % speed of sound
N = conf.N;                    % target length of BRS impulse responses


%% ===== Variables =======================================================
phi = correct_azimuth(phi);


%% ===== BRIR ============================================================
% Initial values
ir_ps = zeros(N,2);

% === Secondary source angle ===
% Calculate the angle between the given loudspeaker and the listener.
% This is needed for the HRIR dataset.
%
%                                 y-axis
%                                    ^
%                                    |
%                X, phi = 0          |
%                 O--------          |
%                  \ a |             |
%                   \ /              |  a = alpha
%                    \               |
%                     o              |
%                     xs             |
%        ----------------------------|--------------------------> x-axis
%
% Angle between listener and source (-pi < alpha <= pi)
% NOTE: phi is the orientation of the listener (see first graph)
[alpha,theta_tmp,r_tmp] = cart2sph(xs(1)-X(1),xs(2)-X(2),0);
% Ensure -pi <= alpha < pi
alpha = correct_azimuth(alpha-phi);

% === IR interpolation ===
ir = get_ir(irs,alpha);
ir_distance = get_ir_distance(irs,alpha);

% === Relative distances and time delay ===
% Time delay of the single source (at the listener position)
% Define an offset to ensure norm(X-xs)-irs.distance+offset > 0
offset = 3; % in m
tau = (norm(X-xs)-ir_distance+offset)/c;
% Time delay in samples
delay = ceil( tau*fs );

% === Amplitude factor ===
% The 1/norm(X-xs) term is for the decreasing of the sound on its way from
% the loudspeaker (xs) to the listener (X). It accounts for the distance that is
% already present in the IR dataset.
amplitude = (1/norm(X-xs)) / (1/ir_distance);

% === Trim IR ===
% make sure the IR has the right length
ir = fix_ir_length(ir,N,delay);

% === Calculate BRIR ===
% Add the point source with the corresponding time delay and amplitude
ir_ps(:,1) = delayline(ir(:,1)',delay,amplitude,conf)';
ir_ps(:,2) = delayline(ir(:,2)',delay,amplitude,conf)';


%% ===== Headphone compensation ==========================================
ir = compensate_headphone(ir_ps,conf);
