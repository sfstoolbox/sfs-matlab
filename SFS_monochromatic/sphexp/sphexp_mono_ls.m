function Anm = sphexp_mono_ls(xs, mode, Nse, f, xq, conf)
%SPHEXP_MONO_LS compute the regular/singular spherical expansion of line source
%
%   Usage: Anm = sphexp_mono_ls(xs, mode, Nse, f, xq, conf)
%
%   Input parameters:
%       xs          - position of line source
%       mode        - 'R' for regular, 'S' for singular
%       Nse         - maximum order of spherical basis functions
%       f           - frequency
%       xq          - optional expansion center coordinate 
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       Al          - regular Spherical Expansion Coefficients
%
%   SPHEXP_MONO_LS(xs, mode, f, Nse, xq, conf) computes the regular/singular 
%   spherical expansion coefficients for a point source at xs. The expansion 
%   will be done around the expansion coordinate xq.
%
%   see also: sphexp_mono_ps sphexp_mono_pw sphexp_convert_circexp

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

%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 6;
narginchk(nargmin,nargmax);
isargposition(xs);
isargpositivescalar(f, Nse);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
    xq = [0,0,0];
else
  isargposition(xq);
end

%% ===== Computation ====================================================
% select suitable basis function
if strcmp('R', mode)
  Am = circexp_mono_ls(xs, mode, Nse, f, xq, conf);
  Anm = sphexp_convert_circexp(Am);
elseif strcmp('S', mode)
  to_be_implemented;
else
  error('%s: %s, unknown mode', upper(mfilename), mode);
end

end

