function irs = fix_irs_length(irs,conf)
%FIX_IRS_LENGTH pads zeros at the beginning of irs set corresponding to its
%maximum claimed distance
%
%   Usage: irs = fix_irs_length(irs,[conf])
%
%   Input parameters:
%       irs  - impulse response data set, e.g. HRTFs
%       conf - optional configuration struct (see SFS_config)
%
%   Output paramteres:
%       irs  - corrected impulse response data set
%
%   FIX_IRS_LENGTH(IRS) pads zeros at the beginning of the given impulse
%   response data set in order to have as many zeros at the beginning as the
%   maximum distance within the data set. This will ensure correct distance
%   extrapolation in get_ir(). Also set the overall length to conf.N if
%   conf.ir.useoriglength is set to "false".
%
%   see also: read_irs, fix_length

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
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    conf = SFS_config;
end
isargstruct(conf);


%% ===== Configuration ===================================================
c = conf.c;
fs = conf.fs;
N = conf.N;


%% ===== Main ============================================================
% get distance of HRTF data set
dist = max(irs.distance);
if dist>10
    warning(['%s: Your maximum distance of the HRTF set is more than 10m. ', ...
        'We will only pad zeros for 10m, this can lead to problems with ', ...
        'get_ir().'],upper(mfilename));
    dist = 10;
end
% append zeros at the beginning of the HRTFs corresponding to its maximum
% distance
samples = ceil(dist/c * fs);
if N-samples<128
    error(['%s: choose a larger conf.N value, because otherwise you ', ...
        'will have only %i samples of your original impulse response.'], ...
        upper(mfilename),N-samples);
end
channels = size(irs.left,2);
irs.left  = [zeros(samples,channels); fix_length(irs.left,N-samples)];
irs.right = [zeros(samples,channels); fix_length(irs.right,N-samples)];
