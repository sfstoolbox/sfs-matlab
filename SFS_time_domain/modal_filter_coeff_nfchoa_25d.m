function [b,a] = modal_filter_coeff_nfchoa_25d(order,R,src,r,conf)
%MODAL_FILTER_COEFF_NFCHOA_25D calculates the modal filter for NFC-HOA 2.5D
%
%   Usage: [b,a] = modal_filter_coeff_nfchoa_25d(order,R,src,[r,[conf]])
%
%   Input parameters:
%       order   - order of filter
%       R       - radius of secondary source distribution (m)
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%       r       - distance of virtual point source (m)
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       b,a     - modal filter coefficients
%
%   MODAL_FILTER_COEFF_NFCHOA_25D(order,R,src,r,conf) returns the modal filter
%   coefficients for 2.5D NFC-HOA synthesis.
%
%   see also: driving_function_imp_nfchoa_25d

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


%% ===== Checking of input  parameters ===================================
nargmin = 3;
nargmax = 5;
narginchk(nargmin,nargmax);
if nargin==nargmax-2
    conf = SFS_config;
elseif nargin==nargmax-1 && ~isstruct(r)
    conf = SFS_config;
end
if strcmp('ps',src) && ~exist('r','var')
    error('%s: the distance r is needed for a point source.',upper(mfilename));
end


%% ===== Configuration ===================================================
c = conf.c;
fs = conf.fs;


%% ===== Computation =====================================================
% compute normalized roots/zeros of spherical Hankel function
% FIXME: this works only until a order of 85. For higher orders factorial will
% return Inf. Hence, for higher orders we have to find another way to calculate
% the zeros of the Hankel funtion.
B=zeros(1,order+2);
A=B;
for n=0:order
        B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
end
B=B(end:-1:1);
% find zeros/roots
z=roots(B);

% compute SOS coefficients of modal driving function
if strcmp('pw',src)
    A(2) = 1;
    p=roots(A);
    [sos,g] = zp2sos(p,z*c/R,2,'down','none');
elseif strcmp('ps',src)
    [sos,g] = zp2sos(z*c/r,z*c/R,1,'up','none');
else
    error('%s: %s is not a known source type.',upper(mfilename),src);
end

% transform coefficients
a = {};
b = {};
for n=1:size(sos,1)
    if isoctave
        [bz,az] = bilinear(sos(n,1:3),sos(n,4:6),1/fs);
    else
        [bz,az] = bilinear(sos(n,1:3),sos(n,4:6),fs);
    end
    b{n} = bz;
    a{n} = az;
end
