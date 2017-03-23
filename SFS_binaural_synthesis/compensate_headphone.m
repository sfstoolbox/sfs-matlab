function ir = compensate_headphone(ir,conf)
%COMPENSATE_HEADPHONE applies a headphone compensation to the IR
%
%   Usage: ir = compensate_headphone(ir,conf)
%
%   Input parameters:
%       ir      - Impulse response to which the compensation should be applied
%       conf    - configuration struct (see SFS_config)
%
%   Output:
%       ir      - Impulse response which is compensated for the given headphone
%
%   COMPENSATE_HEADPHONE(ir,conf) applies a headphone compensation to the
%   given impulse response. Which headphone compensation it should use is
%   mentioned in the conf struct. The compensation filter can be a one-channel
%   (same filter for left and right) or two-channel signal (1st signal: left,
%   2nd signal: right) with signals stored as columns in a matrix.
%   The compensation is only applied, if the conf.ir.usehcomp value is not false.
%
%   See also: ir_wfs, ir_point_source, ir_generic

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
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(ir);
isargstruct(conf);


%% ===== Configuration ===================================================
usehcomp = conf.ir.usehcomp;


%% ===== Computation =====================================================
if(usehcomp)
    lenir = size(ir,1);
    % Read headphone compensation filter
    hcomp = audioread(conf.ir.hcompfile);
    % Check if the IR has the right length for the filter
    if lenir<length(hcomp)
        warning(['The length of the used IR is shorter than the headphone ', ...
            'compensation filter.']);
    end
    % Apply filter
    ir = convolution(hcomp,ir);
    ir = fix_length(ir,lenir);
end
