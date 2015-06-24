function X = sofa_get_listener_position(sofa,coordinate_system)
%SOFA_GET_LISTENER_POSITION returns the listener position from the given SOFA
%data set
%
%   Usage: X = sofa_get_listener_position(sofa,[idx])
%
%   Input parameters:
%       sofa              - impulse response data set (SOFA struct/file)
%       coordinate_system - coordinate system the listener position should be
%                           specified in:
%                             'cartesian' (default)
%                             'spherical'
%
%   Output parameters:
%       X                 - listener position
%
%   SOFA_GET_LISTENER_POSITION(sofa,coordinate_system) returns listener position X
%   as defined in the given SOFA file or struct.
%
%   See also: get_ir, sofa_get_header, sofa_get_listener_view

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
nargmax = 2;
narginchk(nargmin,nargmax)
if nargin==nargmax-1
   coordinate_system = 'cartesian';
else
    isargchar(coordinate_system);
end


%% ===== Computation ====================================================
header = sofa_get_header(sofa);
X = SOFAconvertCoordinates(header.ListenerPosition, ...
                           header.ListenerPosition_Type, ...
                           coordinate_system);
if size(X,1)>1
    error(['%s: you have %i different listener positions, but only one is ', ...
           'currently supported by the SFS Toolbox.'],upper(mfilename), ...
          size(X,1));
end
if strcmp('spherical',coordinate_system)
    X(1) = correct_azimuth(rad(X(1)));
    X(2) = correct_elevation(rad(X(2)));
end
