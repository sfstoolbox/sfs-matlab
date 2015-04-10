function [phi,theta] = sofa_get_head_orientations(sofa,idx)
%SOFA_GET_HEAD_ORIENTATIONS returns phi, theta from the given SOFA data set
%
%   Usage: [phi,theta] = sofa_get_head_orientations(sofa,[idx])
%
%   Input parameters:
%       sofa    - impulse response data set (SOFA struct/file)
%       idx     - index of secondary sources that should be returned.
%                 If no index is specified all sources will be returned.
%
%   Output parameters:
%       phi     - head orientations in the horizontal plane
%       theta   - head orientations in the median plane
%
%   SOFA_GET_HEAD_ORIENTATIONS(sofa,idx) returns head orientation [phi,theta] as
%   defined in the given SOFA file or struct, specified by idx. If no idx is
%   specified all head orientations are returned.
%
%   see also: get_ir, sofa_get_header, sofa_get_secondary_sources

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
