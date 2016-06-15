function D = driving_function_mono_wfs_fs(x0,nx0,xs,f,conf)
%DRIVING_FUNCTION_MONO_WFS_FS returns the driving signal D for a focused source
%in WFS
%
%   Usage: D = driving_function_mono_wfs_fs(x0,nx0,xs,f,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nx0         - directions of the secondary sources / m [nx3]
%       xs          - position of virtual focused source / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_WFS_FS(x0,xs,f,src,conf) returns WFS driving signals
%   for the given secondary sources, the virtual focused source position and the
%   frequency f.
%
%   References:
%       H. Wierstorf, J. Ahrens, F. Winter, F. Schultz, S. Spors (2015) -
%       "Theory of Sound Field Synthesis"
%       S. Spors, H. Wierstorf, M. Geier, J. Ahrens (2009) - "Physical and
%       Perceptual Properties of Focused Sources in Wave Field Synthesis",
%       AES127
%
%   See also: driving_function_mono_wfs, driving_function_imp_wfs_ps

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


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 5;
narginchk(nargmin,nargmax);
isargmatrix(x0,nx0,xs);
isargpositivescalar(f);
isargstruct(conf);
if size(xs,2)~=3
    error('%s: size(xs,2) has to be 3 and not %i.',upper(mfilename),size(xs,2));
end


%% ===== Configuration ==================================================
xref = conf.xref;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% Frequency
omega = 2*pi*f;


if strcmp('2D',dimension) || strcmp('3D',dimension)

    % === 2- or 3-Dimensional ============================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % D using an approximated point sink as source model
        %
        %            1  i w (x0-xs) nx0
        % D(x0,w) = --- --- ------------- e^(i w/c |x0-xs|)
        %           2pi  c  |x0-xs|^(3/2)
        %
        % see Wierstorf et al. (2015), eq.(#D:wfs:fs)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = 1/(2*pi) * 1i*omega/c .* ...
            vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(1i*omega/c.*r);
        %
    elseif strcmp('point_sink',driving_functions)
        % D using a point sink as source model
        %
        % D(x0,w) =
        %
        %  1  /i w      1    \  (x0-xs) nx0
        % --- |--- + ------- |  ----------- e^(i w/c |x0-xs|)
        % 2pi \ c    |x0-xs| /   |x0-xs|^2
        %
        % see Wierstorf et al. (2015), eq.(#D:wfs:fs:woapprox)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = 1/(2*pi) * ( -1i.*omega/c + 1./r ) .* ...
            vector_product(x0-xs,nx0,2) ./ r.^2 .* exp(1i*omega./c.*r);
        %
    elseif strcmp('line_sink',driving_functions)
        % D using a line sink as source model
        %
        %              iw (x0-xs)nk  (1)/ w         \
        % D(x0,w) =  - -- --------- H1  | - |x0-xs| |
        %              2c  |x0-xs|      \ c         /
        %
        % compare Wierstorf et al. (2015), eq.(#D:wfs:fs:ls)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = -1i*omega/(2*c) .* vector_product(x0-xs,nx0,2) ./ r .* ...
            besselh(1,1,omega/c.*r);
        %
    elseif strcmp('delft1988',driving_functions)
        % --- Delft 1988 -------------------------------------------------
        to_be_implemented;
        %
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a focused source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using an approximated point sink as source model
        %                        ___
        %                 g0    |i w (x0-xs) nx0
        % D_2.5D(x0,w) = ---  _ |--- ------------- e^(i w/c |x0-xs|)
        %                2pi   \| c  |x0-xs|^(3/2)
        %
        % see Wierstorf et al. (2015), eq.(#D:wfs:fs:2.5D)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = g0/(2*pi) * sqrt( 1i*omega/c ) .* ...
            vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(1i*omega/c.*r);
        %
    elseif strcmp('spors2009eq7',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using a line sink with point amplitude characteristic as source
        % and the large argument approximation of the driving function above.
        % This results in the "traditional" driving function, derived in
        % Verheijen (1997), see Spors (2009) eq. 7.
        %                       ______
        %                      |   w  |  (x0-xs) nx0
        % D_2.5D(x0,w) = -g0 _ |------  ------------- e^(i w/c |x0-xs|)
        %                     \|2pi ic  |x0-xs|^(3/2)
        %
        % see Spors (2009), eq.(7)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = -g0 * sqrt( omega/(2*pi*c*1i) ) .* ...
            vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(1i*omega/c.*r);
        %
    elseif strcmp('spors2009eq6',driving_functions)
        % --- Spors 2009 -------------------------------------------------
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using a line sink with point source amplitude characteristics as
        % source
        %
        %                    iw (x0-xs) nx0   (1)/ w         \
        % D_2.5D(x0,w) = -g0 -- -----------  H1  | - |x0-xs| |
        %                    2c   |x0-xs|        \ c         /
        %
        % see Spors et al. (2009), eq.(6)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = -g0 * 1i*omega/(2*c) .* ...
            vector_product(x0-xs,nx0,2) ./ r .* besselh(1,1,omega/c.*r);
        %
    elseif strcmp('point_sink',driving_functions)
        % --- Point Sink -------------------------------------------------
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using a point sink as source model
        %
        % D_2.5D(x0,w) =
        %             ___       ___
        %    g0  /   |i w|     | c |    1    \  (x0-xs) nx0
        %   ---  | _ |---  + _ |---  ------- |  ----------- e^(i w/c |x0-xs|)
        %   2pi  \  \| c      \|i w  |x0-xs| /   |x0-xs|^2
        %
        % see Wierstorf et al. (2015), eq.(#D:wfs:fs:woapprox:2.5D)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = g0/(2*pi) .* ( sqrt(1i*omega/c) + sqrt(c/(1i*omega)) ./ r ) .* ...
            vector_product(x0-xs,nx0,2) ./ r.^2 .* exp(1i*omega/c.*r);
        %
    elseif strcmp('delft1988',driving_functions)
        % --- Delft 1988 -------------------------------------------------
        to_be_implemented;
        %
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D focused source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
