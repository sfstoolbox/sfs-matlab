function x0 = sofa_get_emitter_positions(sofa,idx)
%SOFA_GET_EMITTER_POSITIONS returns x0 from the given SOFA data set
%
%   Usage: x0 = sofa_get_emitter_positions(sofa,[idx])
%
%   Input parameters:
%       sofa    - impulse response data set (SOFA struct/file)
%       idx     - index of secondary sources that should be returned.
%                 If no index is specified all sources will be returned.
%
%   Output parameters:
%       x0      - secondary source matrix [n 7]
%
%   SOFA_GET_EMITTER_POSITIONS(sofa,idx) returns secondary sources as defined
%   in the given SOFA file or struct, specified by idx. If no idx is specified
%   all secondary sources are returned.
%
%   see also: sofa_get_header, secondary_source_positions

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

if strcmp('SimpleFreeFieldHRIR',header.GLOBAL_SOFAConventions)
    source_positions = SOFAconvertCoordinates( ...
        header.SourcePosition,header.SourcePosition_Type,'cartesian');
    source_directions = SOFAconvertCoordinates( ...
        header.SourceView,header.SourceView_Type,'cartesian');
    x0 = [source_positions(idx,:) source_directions(idx,:) ...
          ones(size(source_positions(idx,1)))];

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
