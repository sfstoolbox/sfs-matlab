function [d] = driving_function_imp_nfchoa_25d(x0,xs,src,L,conf)
%DRIVING_FUNCTION_IMP_NFCHOA_25D calculates the NFC-HOA 2.5D driving function
%
%   Usage: [d] = driving_function_imp_nfchoa_25d(x0,xs,src,[conf]);
%
%   Input parameters:
%       x0      - position  and direction of secondary sources (m)
%       xs      - position of virtual source or direction of plane wave (m)
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%       L       - diameter of loudspeaker array
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       d  - matrix of driving signals
%
%   DRIVING_FUNCTION_IMP_NFCHOA_25D(x0,xs,src,conf) returns the
%   driving function of 2.5D NFC-HOA for the given source type and position,
%    and loudspeaker positions.
%
%   see also: modal_filter_coeff_nfchoa_25d, wave_field_imp_nfchoa_25d

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


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargsecondarysource(x0)
isargposition(xs);
xs = position_vector(xs);
isargchar(src);
isargpositivescalar(L);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
R = L/2;
nls = size(x0,1);
N = conf.N;
N0 = 0;       % pre-delay

if(isodd(nls))
    order=(nls-1)/2;    %  max order of spherical harmonics
else
    order=floor((nls+1)/2);
end

% if-request as a workaround for the right direction of the wave field
if strcmpi(src,'pw')
    [theta_src, r_src] = cart2pol(-xs(1),xs(2));
elseif strcmpi(src,'ps')
    [theta_src, r_src] = cart2pol(xs(1),-xs(2));
else
    [theta_src, r_src] = cart2pol(xs(1),xs(2));
end

%% ===== Computation =====================================================

% compute impulse responses of modal filters
dm = zeros(order+1,N);
for n=1:order+1
    dm(n,:) = [zeros(1,N0) 1 zeros(1,N-1-N0)];
    [b,a] = modal_filter_coeff_nfchoa_25d(n-1,R,src,r_src,conf);
    for ii=1:length(b)
        dm(n,:) = filter(b{ii},a{ii},dm(n,:));
    end
end

% compute input signal for IFFT
d = zeros(2*order+1,N);
for n=-order:order
    d(n+order+1,:) = dm(abs(n)+1,:) .* exp(1i*n*theta_src);
end

if(iseven(nls))
   d = d(2:end,:);
end

% spatial IFFT
d = circshift(d,[order+1 0]);
d = (2*order+1)*ifft(d,[],1);
d = real(d');
