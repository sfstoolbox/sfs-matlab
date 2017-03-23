function P = norm_sound_field(P,conf)
%NORM_SOUND_FIELD normalizes the sound field
%
%   Usage: P = norm_sound_field(P,conf)
%
%   Input options:
%       P       - sound field
%       conf    - configuration struct (see SFS_config)
%
%   Output options:
%       P       - normalized sound field
%
%   NORM_SOUND_FIELD(P,conf) normalizes the given sound field P. This depends on
%   the conf.plot.normalisation setting. It can be one of the following:
%       'auto'   - if the given absolute sound field value at the center is
%                  > 0.3 it uses automatically 'center', otherwise it uses 'max'
%       'center' - normalises to center of sound field == 1
%       'max'    - normalises to max of sound field == 1
%
%   See also: plot_sound_field

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


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargnumeric(P);
isargstruct(conf);


%% ===== Configuration ===================================================
method = conf.plot.normalisation;


%% ===== Computation =====================================================
if strcmp('auto',method)
    % If sound field at center > 0.3 normalise to center sound field == 1,
    % otherwise normalise to max sound field == 1
    if abs(P(round(end/2),round(end/2)))>0.3
        method = 'center';
    else
        method = 'max';
    end
end
if strcmp('center',method)
    % Center of sound field == 1
    P = P/max(abs(P(round(end/2),round(end/2))));
elseif strcmp('max',method)
    % Max of sound field == 1
    P = P/max(abs(P(:)));
else
    error(['%s: conf.plot.normalisation has to be ''auto'', ''center'' or ', ...
           '''max''.'],upper(mfilename));
end
