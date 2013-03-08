function [d,delay,weight] = driving_function_imp_wfs_25d(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_WFS_25D calculates the WFS 2.5D weighting and delaying
%
%   Usage: [d,delay,weight] = driving_function_imp_wfs_25d(x0,xs,src,[conf]);
%
%   Input parameters:
%       x0      - position  and direction of secondary source (m)
%       xs      - position of virtual source or diirection of plane wave (m)
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%                     'fs' - focused source
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       d       - driving signals
%       delay   - delay of the driving function (s)
%       weight  - weight (amplitude) of the driving function
%
%   DRIVING_FUNCTION_IMP_WFS_25D(x0,xs,src,conf) returns the driving signals and
%   weighting and delay parameters of the WFS 2.5D driving function for the given
%   source type and position and loudspeaker positions.
%
%   see also: wave_field_imp, wave_field_imp_wfs_25d, driving_function_mono_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0)
isargposition(xs);
xs = position_vector(xs);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Speed of sound
c = conf.c;
xref = position_vector(conf.xref);
fs = conf.fs;
usehpre = conf.usehpre;


%% ===== Computation =====================================================

% Calculate pre-equalization filter if required
if usehpre
    hpre = wfs_prefilter(conf);
else
    hpre = 1; % dirac pulse
end

% Secondary source positions and directions
nx0 = x0(:,4:6);
x0 = x0(:,1:3);

% Reference point and source position
xref = repmat(xref,[size(x0,1) 1]);
xs = repmat(xs,[size(x0,1) 1]);

% 2.5D correction factor
g0 = sqrt(2*pi*vector_norm(xref-x0,2));

% Get the delay and weighting factors
if strcmp('pw',src)
    % === Plane wave ===
    % Direction of plane wave
    nxs = bsxfun(@rdivide,xs,vector_norm(xs,2));
    % Delay and amplitude weight
    % NOTE: <n_pw,n(x0)> is the same as the cosinus between their angle
    delay = 1/c * vector_product(nxs,x0,2);
    weight = 2*g0 .* vector_product(nxs,nx0,2);

elseif strcmp('ps',src)
    % === Point source ===
    % Distance between loudspeaker and virtual source
    r = vector_norm(x0-xs,2);
    % Delay and amplitude weight
    delay = r/c;
    weight = g0/(2*pi).*vector_product(x0-xs,nx0,2).*r.^(-3/2);

elseif strcmp('fs',src)
    % === Focused source ===
    % Distance between loudspeaker and virtual source
    r = vector_norm(x0-xs,2);
    % Delay and amplitude weight
    delay =  -r/c;
    weight = g0/(2*pi).*vector_product(x0-xs,nx0,2).*r.^(-3/2);
else
    error('%s: %s is not a known source type.',upper(mfilename),src);
end

% Calculate driving function prototype
% FIXME: check if the zeros at the end are long enough
d_proto = [hpre zeros(1,800)];
d = zeros(length(d_proto),size(x0,1));
for ii=1:size(x0,1)
    % Shift and weight prototype driving function
    % - less delay in driving function is more propagation time in sound
    %   field, hence the sign of the delay has to be reversed in the
    %   argument of the delayline function
    % - the proagation time from the source to the nearest secondary source
    %   is removed
    d(:,ii) = delayline(d_proto,(max(delay)-delay(ii))*fs,weight(ii),conf);
end
