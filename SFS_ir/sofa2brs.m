function brs = sofa2brs(sofa)
%SOFA2BRS converts a sofa file/struct to a brs set suitable for the SSR
%
%   Usage: brs = sofa2brs(sofa)
%
%   Input parameters:
%       sofa    - sofa struct or file name
%
%   Output parameters:
%       brs     - brs data set
%
%   SOFA2BRS(sofa) converts a the given sofa file or struct into a brs set
%   suitable for th SoundScape Renderer. The brs data set is a matrix
%   containing the channels for all directions. As the SoundScape Renderer is
%   currently only working in the horizontal plane, only impulse responses from
%   the sofa data set are used with an elevation of 0.
%
%   See also: sofa_get_header, sofa_get_data, SOFAcalculateAPV

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
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax)


%% ===== Main ===========================================================
header = sofa_get_header(sofa);
if ~strcmp('SimpleFreeFieldHRIR',header.GLOBAL_SOFAConventions)
    error('%s: this SOFA Convention is currently not supported.');
end
% Get available source positions relative to listener view
x0 = SOFAcalculateAPV(header); % / [deg deg m]
% Find elevation = 0
idx = find(abs(x0(:,2))<=eps);
if length(idx)==0
    error(['%s: Your sofa file has no data for an elevation angle of 0deg.' ...
        ' Other angles are not supported by the SoundScape Renderer.'], ...
        upper(mfilename));
end
% Sort azimuth angle in ascending order
[~,sort_idx] = sort(correct_azimuth(rad(x0(idx,1))));
% Get corresponding impulse responses
ir = sofa_get_data(sofa,idx(sort_idx));
% Generate brs set
brs = zeros(sofa.API.N,2*size(ir,1));
for ii = 1:size(ir,1)
    brs(:,ii*2-1:ii*2) = squeeze(ir(ii,:,:))';
end
