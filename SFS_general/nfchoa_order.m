function M = nfchoa_order(nls,conf)
%NFCHOA_ORDER returns the maximum order for the spherical harmonics for the
%given number of secondary sources
%
%   Usage: M = nfchoa_order(nls,conf)
%
%   Input parameters:
%       nls     - number of secondary sources
%
%   Output parameters:
%       M       - spherical harmonics order
%       conf    - configuration struct (see SFS_config)
%
%   NFCHOA_ORDER(nls,conf) returns the maximum order of spherical harmonics for
%   the given number of secondary sources in order to avoid spectral repetitions
%   (spatial aliasing) of the dirving signals. The order is
%
%        / nls/2 - 1,   even nls
%   M = <
%        \ (nls-1)/2    odd nls
%
%   for a circular array and
%         _____
%   M = \|nls/2
%
%   for a spherical array.
%
%   References:
%       J. Ahrens (2012) - "Analytic Methods of Sound Field Synthesis", Springer.
%
%   See also: driving_function_imp_nfchoa, driving_function_mono_nfchoa

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


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargpositivescalar(nls);
isargstruct(conf);


%% ===== Configuration ===================================================
if conf.nfchoa.order
    M = conf.nfchoa.order;
    return;
end
dimension = conf.dimension;


%% ===== Computation =====================================================
% Get maximum order of spherical harmonics to avoid spatial aliasing
if strcmp('2D',dimension) || strcmp('2.5D',dimension)
    % Ahrens (2012), p. 132
    if isodd(nls)
        M = (nls-1)/2;
    else
        M = nls/2 - 1;
    end
elseif strcmp('3D',dimension)
    % Ahrens (2012), p. 125
    M = floor(sqrt(nls/2));
end
