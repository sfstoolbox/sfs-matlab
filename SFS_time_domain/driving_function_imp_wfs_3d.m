function [d,weight,delay] = driving_function_imp_wfs_3d(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_WFS_3D calculates the WFS 3D weighting and delaying
%
%   Usage: [weight,delay] = driving_function_imp_wfs_3d(x0,xs,src,[conf]);
%
%   Input parameters:
%       x0      - position and direction of secondary source (m)
%       xs      - position of virtual source or diirection of plane wave (m)               
%       src     - source type of the virtual source: 'pw' (plane wave), 
%                 'ps' (point source) and 'fs' (focused source)
%                                         
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       weight  - weight (amplitude) of the driving function
%       delay   - delay of the driving function (s)
%       d       - driving function [nx1]
% 
%   DRIVING_FUNCTION_IMP_WFS_3D(x0,xs,src,conf) returns the
%   weighting and delay parameters of the WFS 3D driving function for the given
%   source type and position and loudspeaker positions.
%
%   see also: wave_field_imp_wfs_3d, driving_function_mono_wfs_3d

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

% AUTHOR: Hagen Wierstorf
% $LastChangedDate: 2012-04-24 17:54:07 +0200 (Tue, 24 Apr 2012) $
% $LastChangedRevision: 701 $
% $LastChangedBy: wierstorf.hagen $


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));
% isargsecondarysource(x0)
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
% Check also the activity of the used loudspeaker.

% Calculate pre-equalization filter if required
if usehpre
    hpre = wfs_prefilter3d(conf);
else
    hpre = 1; % dirac pulse
end


% Direction and position of active secondary sources
weights = x0(:,8);
surfaceWeights = x0(:,7);
nx0 = x0(:,4:6);
x0 = x0(:,1:3);

% Reference point and source position
xref = repmat(xref,[size(x0,1) 1]);
xs = repmat(xs,[size(x0,1) 1]);
    
if strcmp('pw',src)
    % === Plane Wave ===
    % Direction of plane wave
    nxs = bsxfun(@rdivide,xs,norm(xs,2));
    % Delay and amplitude weight
    delay = 1/c * vector_product(nxs,x0,2);
    weight = 2 .* vector_product(nxs,nx0,2).* surfaceWeights.* weights;
    
elseif strcmp('ps',src)
    % === Point Source ===
    % Delay and amplitude weight
    delay = vector_norm(x0-xs,2)./c;
    
    weight = (-2.*vector_product(x0-xs,nx0,2)./(vector_norm(x0-xs,2).^2)).*...
             (1./vector_norm(x0-xs,2)+1/c).*surfaceWeights.* weights;

elseif strcmp('fs',src)
% === Focused Source ===  
% Delay and amplitude weight
delay = -1.*vector_norm(x0-xs,2)./c;
weight = (-2.*vector_product(x0-xs,nx0,2)./(vector_norm(x0-xs,2).^2)).*...
         (1./vector_norm(x0-xs,2)+1./c).*surfaceWeights.* weights;

else
    % === Unknown Source Type ===
    error('%s: %s is not a known source type.',upper(mfilename),src);

end

% Calculate driving function prototype
% FIXME: check if the zeros at the end are long enough
delay = delay-min(delay);
d_proto = [hpre zeros(1,800)];
d = zeros(length(d_proto),size(x0,1));
for ii=1:size(x0,1)
    % Shift and weight prototype driving function
    % - less delay in driving function is more propagation time in sound
    %   field, hence the sign of the delay has to be reversed in the
    %   argument of the delayline function
    % - the proagation time from the source to the nearest secondary source
    %   is removed
    d(:,ii) = delayline(d_proto,delay(ii)*fs,weight(ii),conf);
end

end