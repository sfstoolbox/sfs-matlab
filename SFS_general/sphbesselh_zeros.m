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
%   SPHBESSELH_ZEROS(order) finds zeros and roots for a spherical hankel functin
%   of the specified order.
%
%   see also: sphbesselh, driving_function_imp_nfchoa

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


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargpositivescalar(order);


%% ===== Main ============================================================
if order<86
    % --- compute ---
    B = zeros(1,order+2);
    A = B;
    for n=0:order
        B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
    end
    B = B(end:-1:1);
    % find zeros/roots
    z = roots(B);
    % find roots
    A(2) = 1;
    p = roots(A);
else
    % --- use pre-computed ---
    % For orders greater than 85 Matlab/Octave is not able to compute the zeros,
    % because the factorial function returns Inf. We solved this by using the
    % Multiprecission Toolbox from http://www.advanpix.com and the code given at
    % the end of this function. With the Toolbox we were able to compute the
    % zeros up to an order of ... and stored the resulting zeros at
    % http://github.com/sfstoolbox/data/tree/master/sphbesselh_zeros
    filename = sprintf('sphbesselh_zeros_order%04.0f.mat',order);
    file = sprintf('%s/data/sphbesselh_zeros/%s',get_sfs_path(),filename);
    url = ['https://dev.qu.tu-berlin.de/projects/data/repository/revisions/master/' ...
        'raw/sphbesselh_zeros/' filename];
    % download file if not present
    if ~exist(file,'file')
        download_file(url,file);
    end
    load(file);
end
return


%% ===== Computation with Multiprecission Toolbox ========================
B = mp(zeros(1,order+2));
A = B;
for n=mp(0:order)
    B(n+1) = factorial(2*order-n)/(factorial(order-n)*factorial(n)*2^(order-n));
end
B = B(end:-1:1);
z = roots(B);
A(2) = mp(1);
p = roots(A);
