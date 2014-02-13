function [ Lnm ] = asslegendre(n,m,arg)
%ASSLEGENDRE Calculates associates Legendre function of degree n and order m
%
%   Usage: Lnm = asslegendre(n,m,arg)
%
%   Input parameters:
%       n     - spherical harmonic degree
%       m     - spherical harmonic order
%       arg   - values to calculate the lengendre polynomes
%
%   Output parameters:
%       Lnm   - associates Legendre function
%
%   ASSLEGENDRE(n,m,arg) calculates the associates Legendre function of degree
%   n and order m for the values arg.
%
%   see also: sphharmonics, legendre

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
% Copyright (c) 2012      Jens Ahrens                                        *
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


%% ===== Checking input parameters =======================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);
isargpositivescalar(n)
isargscalar(m)
isargnumeric(arg)
if n<abs(m)
    warning( 'Absolute value of order m must be less than or equal to the degree n.' ); 
    Ynm = zeros(size(alpha));
    return;
end


%% ==== Computation ======================================================
Lnm = legendre(n,arg);
if n~=0
    Lnm = squeeze(Lnm(abs(m)+1,:,:));  
end
if m<0
    Lnm = -1.^abs(m) .* factorial(n-abs(m)) ./ factorial(n+abs(m)) .* Lnm;
end
Lnm = reshape(Lnm,size(arg));
