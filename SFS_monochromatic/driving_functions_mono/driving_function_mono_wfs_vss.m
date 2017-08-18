function D = driving_function_mono_wfs_vss(x0,xv,Dv,f,conf)
%DRIVING_FUNCTION_MONO_WFS_VSS returns the driving signal D for a given set of
%virtual secondary sources and their driving signals
%
%   Usage: D = driving_function_mono_wfs_vss(x0,xv,Dv,f,conf)
%
%   Input parameters:
%       x0          - position, direction, and weights of the real secondary
%                     sources / m [nx7]
%       xv          - position, direction, and weights of the virtual secondary
%                     sources / m [mx7]
%       Dv          - driving functions of virtual secondary sources [mx1]
%       f           - frequency of the monochromatic source / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   See also: driving_function_mono_wfs, driving_function_mono_wfs_fs

%   References:
%       S. Spors, J.Ahrens (2010) - "Local Sound Field Synthesis by Virtual
%                                    Secondary Sources", 40th AES

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
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
isargvector(Dv);
isargpositivescalar(f);
isargsecondarysource(x0,xv);
isargstruct(conf);


%% ===== Configuration ==================================================
dimension = conf.dimension;


%% ===== Computation ====================================================
% Get driving signals
if strcmp('2.5D',dimension) || strcmp('3D',dimension)
    % === Focussed Point Sink ===
    conf.driving_functions = 'default';
elseif strcmp('2D',dimension)
    % === Focussed Line Sink ===
    % We have to use the driving function setting directly, because in opposite
    % to the case of a non-focused source where 'ps' and 'ls' are available as
    % source types, for a focused source only 'fs' is available.
    % Have a look at driving_function_mono_wfs_fs() for details on the
    % implemented focused source types.
    conf.driving_functions = 'line_sink';
else
    error('%s: %s is not a known source type.',upper(mfilename),dimension);
end

% Get driving signals for real secondary sources
%
% See Spors (2010), fig. 2 & eq. (12)
Ns = size(xv,1);
N0 = size(x0,1);

Dmatrix = zeros(N0,Ns);

for idx=1:Ns
    [xtmp, xdx] = secondary_source_selection(x0,xv(idx,1:6),'fs');
    if (~isempty(xtmp))
        wtap = tapering_window(xtmp,conf);
        Dmatrix(xdx,idx) = ...
            driving_function_mono_wfs(xtmp,xv(idx,1:3),'fs',f,conf) .* wtap;
    end
end

D = Dmatrix*(Dv.*xv(:,7));
