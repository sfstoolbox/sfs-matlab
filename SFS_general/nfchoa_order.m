function M = nfchoa_order(nls,conf)
%NFCHOA_ORDER returns the maximum order for the spherical harmonics for the
%given number of secondary sources
%
%   Usage: M = nfchoa_order(nls,[conf])
%
%   Input parameters:
%       nls     - number of secondary sources
%
%   Output parameters:
%       M       - spherical harmonics order
%       conf    - optional configuration struct (see SFS_config)
%
%   NFCHOA_ORDER(nls) returns the maximum order of spherical harmonics for the
%   given number of secondary sources in order to avoid spectral repetitions
%   (spatial aliasing) of the dirving signals. The order is
%
%        / nls/2,       even nls
%   M = <
%        \ (nls-1)/2    odd nls
%
%   for a circular array and 
%         _____
%   M = \|nls/2
%
%   for a spherical array.
%
%   References:
%       J. Ahrens (2012) - "Analytic Methods of Sound Field Synthesis", Springer.
%
%   see also: driving_function_imp_nfchoa, driving_function_mono_nfchoa

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
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    conf = SFS_config;
end
isargpositivescalar(nls);
isargstruct(conf);


%% ===== Configuration ===================================================
if conf.nfchoa.order
    M = conf.nfchoa.order;
    return;
end
dimension = conf.dimension;


%% ===== Computation =====================================================
% get maximum order of spherical harmonics to avoid spatial aliasing
if strcmp('2D',dimension) || strcmp('2.5D',dimension)
    % Ahrens (2012), p. 132
    if isodd(nls)
        M = (nls-1)/2;
    else
        M = nls/2;
    end
elseif strcmp('3D',dimension)
    % Ahrens (2012), p. 125
    M = floor(sqrt(nls/2));
end
