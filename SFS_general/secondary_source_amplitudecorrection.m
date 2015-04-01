function x0 = secondary_source_amplitudecorrection(x0)
%SECONDARY_SOURCE_AMPLITUDECORRECTION the secondary sources
%
%   Usage: x0 = secondary_source_amplitudecorrection(x0)
%
%   Input parameters:
%       x0          - secondary sources / m
%
%   Output options:
%       x0          - secondary sources / m, containing the applied amplitude 
%                     correction in its weights in x0(:,7)
%
%   SECONDARY_SOURCE_AMPLITUDECORRECTION(x0)

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax),
isargsecondarysource(x0);

%% ===== Calculation ====================================================
if size(x0,1) < 3
  % nothing to normalize for this
  return;
end

xd = vector_norm(x0(2:end,1:3) - x0(1:end-1,1:3),2);
xd = xd./sum(xd,1);

x0(1,7) = xd(1);
x0(2:end-1,7) = xd(1:end-1) + xd(2:end);
x0(end,7) = xd(end);

x0(:,7) = x0(:,7)./2;

end