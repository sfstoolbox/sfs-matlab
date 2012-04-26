function [hd] = modal_filter_pw_nfchoa_25d(R,order,fs)
%MODAL_FILTER_PW_NFCHOA_25D calculates the modal filter for NFC-HOA 2.5D
%
%   Usage: [hd] = modal_filter_pw_nfchoa_25d(R,order,fs)
%
%   Input parameters:
%       R       - radius of secondary source distribution (m)
%       order   - order of filter
%       fs      - sampling frequency
%
%   Output parameters:
%       hd  - modal filter object
%
%   MODAL_FILTER_PW_NFCHOA_25D(R,order,fs) returns the modal filter
%   object for 2.5D NFC-HOA synthesis of a plane wave.
%
%   see also: modal_filter_ps_nfchoa_25d

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

% AUTHOR: Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


c=343;

% compute normalized roots/zeros of spherical Hankel function
B=zeros(1,order+2);
A=B;

for n=0:order
        B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
end
B=B(end:-1:1);
%A(2) = (-1)^order;
A(2) = 1;

% find zeros/roots
z=roots(B);
p=roots(A);


% compute SOS coefficients of modal driving function
[sos,g] = zp2sos(p,z*c/R,2,'down', 'none');

% transform coefficients
for n=1:size(sos,1)
    %[bz,az] = impinvar(sos(n,1:3),sos(n,4:6),fs);
    [bz,az] = bilinear(sos(n,1:3),sos(n,4:6),fs);
    if(length(bz)==2)
        sos(n,2:3)=bz;
        sos(n,5:6)=az;
    else
        sos(n,1:3)=bz;
        sos(n,4:6)=az;
    end
end

if isoctave
    error(['%s: the HOA implementation depends on dfilt at the moment, ', ...
        'which is not available under Octave'],upper(mfilename));
else
    % realize FOS/SOS as DF-II structure
    hd = dfilt.df2sos(sos); 
    hd.ScaleValues(end)=2*(-1)^order;
end
