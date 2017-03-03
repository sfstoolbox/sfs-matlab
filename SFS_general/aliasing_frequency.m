function [fal,dx0] = aliasing_frequency(x0,conf)
%ALIASING_FREQUENCY returns the aliasing frequency for the given secondary
%sources
%
%   Usage: [fal,dx0] = aliasing_frequency([x0],conf)
%
%   Input options:
%       x0      - secondary sources / m
%       conf    - configuration struct (see SFS_config)
%
%   Output options:
%       fal     - aliasing frequency / Hz
%       dx0     - mean distance between secondary sources / m
%
%   ALIASING_FREQUENCY(x0,conf) returns the aliasing frequency for the given
%   secondary sources. First the mean distance dx0 between the secondary sources
%   is calculated, afterwards the aliasing frequency is calculated after Spors
%   (2009) as fal = c/(2*dx0). If no secondary sources x0 are provided, they are
%   first calculated by calling secondary_source_positions().
%   For a calculation that includes the dependency on the listener position have
%   a look at Start (1997).
%
%   References:
%       S. Spors and J. Ahrens (2009) - "Spatial sampling artifacts of wave field
%       synthesis for the reproduction of virtual point sources", 126th AES Conv.
%       E. Start (1997) - "Direct Sound Enhancement by Wave Field Synthesis",
%       TU Delft.
%
%   See also: sound_field_mono_wfs, secondary_source_positions,
%       secondary_source_distance

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
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = x0;
    x0 = [];
end
isargstruct(conf);


%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation =====================================================
% If no explicit secondary source distribution is given, calculate one
if isempty(x0)
    x0 = secondary_source_positions(conf);
end
% Get average distance between secondary sources
dx0 = secondary_source_distance(x0);
% Calculate aliasing frequency
fal = c/(2*dx0);
