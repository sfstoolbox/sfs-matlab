function [P, x, y, z] = sound_field_mono_sphexp(X,Y,Z,ABnm,mode,f,xq,conf)
%SOUND_FIELD_MONO_SPHEXPR simulates a sound field with regular spherical
%expansion coefficients
%
%   Usage: [P, x, y, z] = sound_field_mono_sphexp(X,Y,Z,Al,f,x0,conf)
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax]
%       Y           - y-axis / m; single value or [ymin,ymax]
%       Z           - z-axis / m; single value or [zmin,zmax]
%       ABnm        - regular/singular spherical expansion coefficients
%       mode        - 'R' for regular, 'S' for singular
%       f           - frequency in Hz
%       xq          - optional expansion center coordinates, default: [0, 0, 0]
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - resulting soundfield
%
%   SOUND_FIELD_MONO_SPHEXPR(X,Y,Z,ABnm,mode,f,xq,conf)
%
%   see also: sphbasis_mono_XYZgrid, sound_field_mono_basis

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
nargmin = 6;
nargmax = 8;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z,ABnm);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
  xq = [0, 0, 0];
end  
isargposition(xq);

%% ===== Computation ====================================================
[jn, hn, Ynm, x, y, z] = sphbasis_mono_XYZgrid(X,Y,Z,f,xq,conf);

if strcmp('R', mode)
  P = sound_field_mono_basis(ABnm,jn,Ynm,conf);
elseif strcmp('S', mode)
  P = sound_field_mono_basis(ABnm,hn,Ynm,conf);
else
  error('unknown mode:');
end
  
end

