function X = sofa_get_listener_position(sofa,coordinate_system)
%SOFA_GET_LISTENER_POSITION returns the listener position from the given SOFA
%data set
%
%   Usage: X = sofa_get_listener_position(sofa,[coordinate_system])
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
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
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
