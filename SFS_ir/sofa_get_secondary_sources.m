function [x0,nss] = sofa_get_secondary_sources(sofa,idx,coordinate_system)
%SOFA_GET_SECONDARY_SOURCES returns x0 from the given SOFA data set
%
%   Usage: [x0,nss] = sofa_get_secondary_sources(sofa,[idx],[coordinate_system])
%
%   Input parameters:
%       sofa              - impulse response data set (SOFA struct/file)
%       idx               - index of secondary sources that should be returned
%                           (default: return all secondary sources).
%       coordinate_system - coordinate system the position and direction of the
%                           secondary sources should be specified in:
%                             'cartesian' (default)
%                             'spherical'
%
%   Output parameters:
%       x0      - secondary source matrix [n 7]
%       nss     - number of secondary sources
%
%   SOFA_GET_SECONDARY_SOURCES(sofa,idx,coordinate_system) returns secondary
%   sources as defined in the given SOFA file or struct, specified by idx.
%   If no idx is specified all secondary sources are returned. The coordinate
%   system the position and direction of the secondary sources are specified in
%   can be given by the string 'coordinate_system', 'cartesian' is assumed as
%   default.
%
%   see also: sofa_get_header, secondary_source_positions

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
nargmax = 3;
narginchk(nargmin,nargmax)
if nargin==nargmax-2
    idx = ':';
    coordinate_system = 'cartesian';
elseif nargin==nargmax-1
    if ischar(idx)
        coordinate_system = idx;
        idx = ':';
    else
        coordinate_system = 'cartesian';
    end
else
    isargvector(idx);
    isargchar(coordinate_system);
end


%% ===== Computation ====================================================
header = sofa_get_header(sofa);

if strcmp('SimpleFreeFieldHRIR',header.GLOBAL_SOFAConventions)
    % For free field HRTFs the source positions are equivalent to the apparent
    % positons of the sources
    apparent_directions = SOFAcalculateAPV(header);
    [apparent_positions(:,1) apparent_positions(:,2) apparent_positions(:,3)] = ...
        sph2cart(rad(apparent_directions(:,1)), ...
                 rad(apparent_directions(:,2)), ...
                 apparent_directions(:,3));
    x0 = [apparent_positions(idx,:) ...
          direction_vector(apparent_positions(idx,:), ...
                           zeros(size(apparent_positions,1),3)) ...
          ones(size(apparent_positions(idx,1)))];

elseif strcmp('SingleRoomDRIR',header.GLOBAL_SOFAConventions)
    to_be_implemented;

elseif strcmp('MultiSpeakerBRIR',header.GLOBAL_SOFAConventions)
    source_position = SOFAconvertCoordinates( ...
        header.SourcePosition,header.SourcePosition_Type,'cartesian');
    emitter_positions = SOFAconvertCoordinates( ...
        header.EmitterPosition,header.EmitterPosition_Type,'cartesian');
    emitter_directions = SOFAconvertCoordinates( ...
        header.EmitterView,header.EmitterView_Type,'cartesian');
    x0 = [bsxfun(@minus,emitter_positions(idx,:),source_position) ...
          emitter_directions(idx,:) ...
          ones(size(emitter_positions(idx,1)))];

else
    error('%s: %s convention currently not supported.', ...
        upper(mfilename),header.GLOBAL_SOFAConventions);
end

nss = size(x0,1);

if strcmp('cartesian',coordinate_system)
    return;
elseif strcmp('spherical',coordinate_system)
    [x0(:,1) x0(:,2) x0(:,3)] = cart2sph(x0(:,1),x0(:,2),x0(:,3));
else
    error('%s: %s is not a supported coordinate system.', ...
        upper(mfilename),coordinate_system);
end
