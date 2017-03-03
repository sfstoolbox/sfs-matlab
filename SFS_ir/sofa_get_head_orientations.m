function [phi,theta,r] = sofa_get_head_orientations(sofa,idx)
%SOFA_GET_HEAD_ORIENTATIONS returns phi, theta and r from the given SOFA data set
%
%   Usage: [phi,theta,r] = sofa_get_head_orientations(sofa,[idx])
%
%   Input parameters:
%       sofa   - impulse response data set (SOFA struct/file)
%       idx    - index of measurement for which head orientation should be
%                returned (default: return all head orientations)
%
%   Output parameters:
%       phi    - head orientations in the horizontal plane / rad
%       theta  - head orientations in the median plane / rad
%       r      - head orientations radii / m
%
%   SOFA_GET_HEAD_ORIENTATIONS(sofa,idx) returns head orientation [phi,theta,r] as
%   defined in the given SOFA file or struct, specified by idx. If no idx is
%   specified, all head orientations are returned.
%
%   See also: get_ir, sofa_get_header, sofa_get_secondary_sources

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
    idx = ':';
else
    isargvector(idx);
end


%% ===== Computation ====================================================
header = sofa_get_header(sofa);
listener_view = SOFAconvertCoordinates(header.ListenerView, ...
                                       header.ListenerView_Type,'spherical');
phi = correct_azimuth(rad(listener_view(idx,1)));
theta = correct_elevation(rad(listener_view(idx,2)));
r = listener_view(idx,3);
