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
%       conf        - configuration struct (see SFS_config)
%
%   Output paramteres:
%       ir          - impulse responses [M C N]
%
%   IR_CORRECT_DISTANCE(ir,ir_distance,r,conf) weights and delays the given
%   impulse responses, that they are conform with the specified distance r.
%   The impulse responses are zero-padded to not loose information.
%
%   See also: get_ir

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
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


%% ===== Configuration ==================================================
c = conf.c;
fs = conf.fs;


%% ===== Computation ====================================================
% Time delay to ensure zeros before the maximum radius
offset = ceil(ir_distance/c*fs); % / samples
% Time delay of the source (at the listener position)
delay = (r-ir_distance)/c*fs; % / samples
% Amplitude weighting (point source model)
% This gives weight=1 for r==ir_distance
weight = ir_distance./r;
% Make sure impulse response is long enough for the specified delay
ir = cat(3,ir,zeros(size(ir,1),size(ir,2),max(ceil(delay+offset))));
% Apply delay and weighting
ir = delayline(ir,[delay+offset; delay+offset],[weight; weight],conf);
