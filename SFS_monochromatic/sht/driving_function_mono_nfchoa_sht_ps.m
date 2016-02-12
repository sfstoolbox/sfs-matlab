function Dnm = driving_function_mono_nfchoa_sht_ps(xs, Nse, f, conf)
%DRIVING_FUNCTION_MONO_NFCHOA_SHT_PS computes the spherical harmonics 
%transform of nfchoa driving functions for a point source.
%
%   Usage: D = driving_function_mono_nfchoa_sht_ps(xs, Nse, f, conf)
%
%   Input parameters:
%       xs          - position of point source [1 x 3] / m
%       Nse         - maximum order of spherical basis functions
%       f           - frequency [m x 1] or [1 x m] / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       Dnm         - regular spherical harmonics transform of driving
%                     function signal [n x m]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_SHT_PS(xs, Nse, f, conf) returns spherical 
%   harmonics transform of the NFCHOA driving function with maximum order Nse
%   for a virtual point source at xs.
%
%   see also: sphexp_mono_ps driving_function_mono_nfchoa_sht_sphexp

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

%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargposition(xs);
isargpositivescalar(Nse);
isargvector(f);
isargstruct(conf);

%% ===== Configuration ==================================================
X0 = conf.secondary_sources.center;

%% ===== Computation ====================================================
% Calculate spherical expansion coefficients of point source
Pnm = sphexp_mono_ps(xs, 'R', Nse, f, X0, conf);
% Calculate spherical harmonics driving function
Dnm = driving_function_mono_nfchoa_sht_sphexp(Pnm, f, conf);
