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
