function D = driving_function_mono_wfs_pw(x0,nx0,nk,f,conf)
%DRIVING_FUNCTION_MONO_WFS_PW driving signal for a plane wave in WFS
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
%   See also: driving_function_mono_wfs, driving_function_imp_wfs_ps
%
%   References:
%       Schultz (2016) - "Sound Field Synthesis for Line Source Array
%       Applications in Large-Scale Sound Reinforcement", PhD thesis,
%       Universität Rostock, https://doi.org/10.18453/rosdok_id00001765

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
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
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
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

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % D using a plane wave as source model
        %
        %              i w
        % D(x0,w) =  2 --- nk nx0  e^(-i w/c nk x0)
        %               c
        %
        % https://sfs.rtfd.io/en/3.2/d_wfs/#equation-fd-wfs-plane
        %
        D = 2.*1i.*omega./c .* vector_product(nk,nx0,2) ...
            .* exp(-1i.*omega./c.*vector_product(nk,x0,2));
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);

    switch driving_functions
    case {'default', 'reference_point'}
        % Driving function with only one stationary phase approximation, i.e.
        % reference to one point in field
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt( 2*pi*vector_norm(xref-x0,2) );
        %
        % D_2.5D using a plane wave as source model
        %                              ___
        %                             |i w|
        % D_2.5D(x0,w) = 2g0 nk nx0 _ |---  e^(-i w/c nk x0)
        %                            \| c
        %
        % https://sfs.rtfd.io/en/3.2/d_wfs/#equation-fd-wfs-plane-25d
        %
        D = 2.*g0 .* vector_product(nk,nx0,2) .* sqrt(1i.*omega./c) ...
            .* exp(-1i.*omega./c.*vector_product(nk,x0,2));
        %
    case {'reference_line'}
        % Driving function with two stationary phase approximations,
        % reference to a line parallel to a LINEAR secondary source distribution
        %
        % Distance ref-line to linear ssd
        dref = vector_product(xref-x0,nx0,2);
        %
        % 2.5D correction factor
        %        ____________
        % g0 = \| 2pi * d_ref
        %
        g0 = sqrt(2.*pi.*dref);
        %
        % D_2.5D using a plane wave as source model
        %                                ___
        %                      ______   |i w|
        % D_2.5D(x0,w) = 2g0 \|nk nx0 _ |---  e^(-i w/c nk x0)
        %                              \| c
        %
        % See Schultz (2016), eq. (2.170)
        %
        D = 2.*g0 .* sqrt(1i*omega/c) .* sqrt(vector_product(nk,nx0,2)) ...
            .* exp(-1i.*omega./c.*vector_product(nk,x0,2));
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
          'for a 2.5D plane wave.'],upper(mfilename),driving_functions);
    end
else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
