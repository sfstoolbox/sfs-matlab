function [z,p] = sphbesselh_zeros(order)
%SPHBESSELH_ZEROS finds zeros/roots of spherical hankel function
%
%   Usage: [z,p] = sphbesselh_zeros(order)
%
%   Input parameters:
%       order       - order of hankel function
%
%   Output parameters:
%       z       - zeros/roots fo the Bessel function
%       p       - roots of the Bessel function
%
%   SPHBESSELH_ZEROS(order) finds zeros and roots for a spherical hankel function
%   of the specified order. Due to numerical problems, the order is limited up
%   to 85.
%
%   See also: sphbesselh, driving_function_imp_nfchoa

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargpositivescalar(order);


%% ===== Main ============================================================
if order<86
    % Formula for nominator (source?)
    B = zeros(1,order+2);
    for n=0:order
        B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
    end
    B = B(end:-1:1);
    % Zeros
    z = roots(B);
    % Poles (are always zero)
    p = zeros(order,1);
else
    error(['%s: for orders higher than 85 no stable numerical ', ...
           'method is available at the moment to caclulate the zeros.'], ...
          upper(mfilename));
end
return


%% ===== Computation with Multiprecission Toolbox ========================
% For the Multiprecission Toolbox, see: http://www.advanpix.com
% Unfortunately it turned out, that the obtained zeros with this method have
% some systematic errors, see
% https://github.com/sfstoolbox/sfs/issues/57#issuecomment-183791477
% The following code was used to calculate the zeros with the Multiprecission
% Toolbox. The results are stored at
% https://github.com/sfstoolbox/data/tree/master/sphbesselh_zeros
B = mp(zeros(1,order+2));
A = B;
for n=mp(0:order)
    B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
end
B = B(end:-1:1);
z = roots(B);
A(2) = mp(1);
p = roots(A);
