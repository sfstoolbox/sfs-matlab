function ir = ir_correct_distance(ir,ir_distance,r,conf)
%IR_CORRECT_DISTANCE weights and delays the impulse reponse for a desired
%distance
%
%   Usage: ir = ir_correct_distance(ir,ir_distance,r,conf)
%
%   Input parameters:
%       ir          - impulse responses [M C N]
%                       M ... number of measurements
%                       C ... number of receiver channels
%                       N ... number of samples
%       ir_distance - distance of the given impulse responses [M 1]
%       r           - desired distance [1]
%       conf        - configuration struct (see SFS_config), containing:
%                       conf.c
%                       conf.fs
%                       conf.N
%                       conf.ir.useoriglength
%
%   Output paramteres:
%       ir          - impulse responses [M C N]
%
%   IR_CORRECT_DISTANCE(ir,ir_distance,r,conf) weights and delays the given
%   impulse responses, that they are conform with the specified distance r.
%
%   See also: get_ir

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
%nargmin = 4;
%nargmax = 4;
%narginchk(nargmin,nargmax);


%% ===== Configuration ==================================================
c = conf.c;
fs = conf.fs;
useoriglength = conf.ir.useoriglength;
%N = conf.N; is used if useoriglength==false


%% ===== Computation ====================================================
% Stop extrapolation for distances larger than 10m
if any(ir_distance>10)
    ir_distance = min(ir_distance,10);
    warning(['%s: Your desired radius is larger than 10m, but we will ', ...
        'only extrapolate up to 10m. All larger radii will be set to ', ...
        '10m.'],upper(mfilename));
end
% Append zeros at the beginning of the impulse responses corresponding to
% its maximum radius
if ~useoriglength
    zero_padding = ceil(ir_distance/c * fs);
    if conf.N-zero_padding<128
        error(['%s: choose a larger conf.N value, because otherwise you ', ...
            'will have only %i samples of your original impulse response.'], ...
            upper(mfilename),conf.N-zero_padding);
    end
else
    zero_padding = 0;
end
% Time delay of the source (at the listener position)
delay = (r-ir_distance)/c*fs; % / samples
% Amplitude weighting (point source model)
% This gives weight=1 for r==ir_distance
weight = ir_distance./r;
if abs(delay)>size(ir,3)
    error(['%s: your impulse response is to short for a desired ', ...
        'delay of %.1f samples.'],upper(mfilename),delay);
end
% Apply delay and weighting
ir = delayline(ir,[delay+zero_padding; delay+zero_padding], ...
               [weight; weight],conf);
