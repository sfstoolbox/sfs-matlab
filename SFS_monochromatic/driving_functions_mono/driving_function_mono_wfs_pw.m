function D = driving_function_mono_wfs_pw(x0,nx0,nk,f,conf)
%DRIVING_FUNCTION_MONO_WFS_PW returns the driving signal D for a plane wave in
%WFS
%
%   Usage: D = driving_function_mono_wfs_pw(x0,nx0,nk,f,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nx0         - directions of the secondary sources / m [nx3]
%       nk          - direction of plane wave / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_WFS_PW(x0,nx0,nk,f,conf) returns WFS driving signals
%   for the given secondary sources, the virtual plane wave direction and the
%   frequency f.
%
%   References:
%       H. Wierstorf, J. Ahrens, F. Winter, F. Schultz, S. Spors (2015) -
%       "Theory of Sound Field Synthesis"
%
%   See also: driving_function_mono_wfs, driving_function_imp_wfs_ps

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Team                                   *
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
isargmatrix(x0,nx0,nk);
isargpositivescalar(f);
isargstruct(conf);


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
        % D using a plane wave as source model
        %
        %              i w
        % D(x0,w) =  2 --- nk nx0  e^(-i w/c nk x0)
        %               c
        %
        % see Wierstorf et al. (2015), eq.(#D:wfs:pw)
        %
        D = 2*1i*omega/c .* vector_product(nk,nx0,2) .*  ...
            exp(-1i*omega/c.*vector_product(nk,x0,2));
        %
    elseif strcmp('delft1988',driving_functions)
        % --- Delft 1988 -------------------------------------------------
        to_be_implemented;
        %
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a plane wave.'],upper(mfilename),driving_functions);
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
        % D_2.5D using a plane wave as source model
        %                               ___
        %                              | w |
        % D_2.5D(x0,w) = 2 g0 nk nx0 _ |---  e^(-i w/c nk x0)
        %                             \|i c
        %
        % see Wierstorf et al. (2015), eq.(#D:wfs:pw:2.5D)
        %
        D = 2*g0 .* vector_product(nk,nx0,2) .* sqrt(omega/(1i*c)) .* ...
            exp(-1i*omega/c.*vector_product(nk,x0,2));
        %
    elseif strcmp('delft1988',driving_functions)
        % --- Delft 1988 -------------------------------------------------
        to_be_implemented;
        %
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D plane wave.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
