function D = driving_function_mono_wfs_ps(x0,nx0,xs,f,conf)
%DRIVING_FUNCTION_MONO_WFS_PS driving signal for a point source in WFS
%
%   Usage: D = driving_function_mono_wfs_ps(x0,nx0,xs,f,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nx0         - directions of the secondary sources / m [nx3]
%       xs          - position of virtual point source / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   See also: driving_function_mono_wfs, driving_function_imp_wfs_ps
%
%   References:
%       Spors and Ahrens (2010) - "Analysis and Improvement of Pre-equalization
%       in 2.5-Dimensional Wave Field Synthesis", 128th Convention of the Audio
%       Engineering Society, Paper 8121,
%       http://www.aes.org/e-lib/browse.cfm?elib=15418
%
%       Spors, Rabenstein, Ahrens (2008) - "The Theory of Wave Field Synthesis
%       Revisited", 124th Convention of the Audio Engineering Society, Paper
%       7358, http://www.aes.org/e-lib/browse.cfm?elib=14488
%
%       Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, TU Berlin,  https://doi.org/10.14279/depositonce-4310

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
isargmatrix(x0,nx0,xs);
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
        % D using a point source and the approximation w/c|x0-xs|>>1
        %
        %            1  i w (x0-xs) nx0
        % D(x0,w) = --- --- ----------- e^(-i w/c |x0-xs|)
        %           2pi  c   |x0-xs|^2
        %
        % https://sfs.rtfd.io/en/3.2/d_wfs/#equation-fd-wfs-point
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = 1./(2.*pi) .* (1i.*omega)./c ...
            .* vector_product(x0-xs,nx0,2) ./ r.^2 ...
            .* exp(-1i.*omega./c.*r);
        %
    case 'point_source'
        % D using a point source as source model
        %
        %            1  / i w      1    \  (x0-xs) nx0
        % D(x0,w) = --- | --- - ------- |  ----------- e^(-i w/c |x0-xs|)
        %           2pi \  c    |x0-xs| /   |x0-xs|^2
        %
        % https://sfs.rtfd.io/en/3.2/d_wfs/#equation-fd-wfs-point-woapprox
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = 1./(2.*pi) .* ( (1i.*omega)./c - 1./r ) ...
            .* vector_product(x0-xs,nx0,2) ./ r.^2 ...
            .* exp(-1i.*omega./c.*r);
        %
    case 'line_source'
        % D using a line source as source model (truly 2D model)
        %
        %             1  i w (x0-xs) nx0  (2)/ w         \
        % D(x0,x) = - -- --- ----------- H1  | - |x0-xs| |
        %             2c  c    |x0-xs|       \ c         /
        %
        % See Spors et al. (2008), eq.(23)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = -1./(2.*c) .* 1i.*omega./c ...
            .* vector_product(x0-xs,nx0,2) ./ r ...
            .* besselh(1,2,omega./c.*r);
        %
    case 'legacy'
        % --- Old SFS Toolbox default ------------------------------------
        % D using a point sink and large distance approximation
        %
        %            1  i w  (x0-xs) nx0
        % D(x0,w) = --- --- ------------- e^(-i w/c |x0-xs|)
        %           2pi  c  |x0-xs|^(3/2)
        %
        % See Wierstorf (2014), eq.(2.61)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = 1./(2.*pi) .* (1i.*omega)./c ...
            .* vector_product(x0-xs,nx0,2) ./ r.^(3./2) ...
            .* exp(-1i.*omega./c.*r);
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);

    switch driving_functions
    case {'default', 'reference_point', 'opperschall', 'volk2010'}
        % Driving function with only one stationary phase approximation,
        % reference to one point in field
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % 2.5D correction factor
        %         _____________________
        %        |      |xref-x0|
        % g0 = _ |---------------------
        %       \| |xref-x0| + |x0-xs|
        %
        g0 = sqrt( vector_norm(xref-x0,2) ./ (vector_norm(xref-x0,2) + r) );
        %
        % D_2.5D(x0,w) =
        %       ___    ___
        %      | 1    |i w (x0-xs) nx0
        % g0 _ |--- _ |--- ------------- e^(-i w/c |x0-xs|)
        %     \|2pi  \| c  |x0-xs|^(3/2)
        %
        % https://sfs.rtfd.io/en/3.2/d_wfs/#equation-fd-wfs-point-25d
        %
        % Driving signal
        D = 1./sqrt(2*pi) .* sqrt(1i.*omega./c) .* g0 ...
          .* vector_product(x0-xs,nx0,2) ./ r.^(3./2) ...
          .* exp(-1i.*omega./c.*r);
        %
    case {'reference_line', 'delft1988'}
        % Driving function with two stationary phase approximations,
        % reference to a line parallel to a LINEAR secondary source distribution
        %
        % Distance ref-line to linear ssd
        dref = abs( vector_product(xref-x0,nx0,2) );
        % Distance source and linear ssd
        ds = abs( vector_product(xs-x0,nx0,2) );
        %
        % 2.5D correction factor
        %        _______________________
        % g0 = \| d_ref / (d_ref + d_s)
        %
        g0 = sqrt( dref ./ (dref + ds) );
        %
        % D_2.5D(x0,w) =
        %       ___    ___
        %      | 1    |i w (x0-xs) nx0
        % g0 _ |--- _ |--- ------------- e^(-i w/c |x0-xs|)
        %     \|2pi  \| c  |x0-xs|^(3/2)
        %
        % https://sfs.rtfd.io/en/3.2/d_wfs/#equation-fd-wfs-point-25d-refline
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = 1./sqrt(2*pi) .* sqrt(1i.*omega./c) .* g0 ...
          .* vector_product(x0-xs,nx0,2) ./ r.^(3./2) ...
          .* exp(-1i.*omega./c.*r);
        %
    case 'legacy'
        % --- Old SFS Toolbox default ------------------------------------
        % 2.5D correction factor
        %        _____________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using a point source and large distance approximation
        %                       ___
        %                g0    |i w  (x0-xs) nx0
        % D_2.5D(x0,w) = --- _ |--- ------------- e^(-i w/c |x0-xs|)
        %                2pi  \| c  |x0-xs|^(3/2)
        %
        % See Wierstorf (2014), eq.(2.62)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = g0./(2.*pi) .* sqrt(1i.*omega./c) ...
            .* vector_product(x0-xs,nx0,2) ./ r.^(3./2) ...
            .* exp(-1i.*omega./c.*r);
        %
    case 'point_source'
        % --- Wierstorf 2014, without approximation ----------------------
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using a point source as source model
        %
        % D_2.5D(x0,w) =
        %           ___       ___
        %  g0  /   |i w      | c      1    \  (x0-xs) nx0
        % ---  | _ |---  - _ |---  ------- |  ----------- e^(-i w/c |x0-xs|)
        % 2pi  \  \| c      \|i w  |x0-xs| /   |x0-xs|^2
        %
        % See Wierstorf (2014), eq.(2.60)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = g0./(2.*pi) .* ( sqrt(1i.*omega./c) - sqrt(c./(1i.*omega) ./ r ) ) ...
            .* vector_product(x0-xs,nx0,2) ./ r.^2 ...
            .* exp(-1i.*omega./c .* r);
        %
    case 'SDMapprox'
        % --- Spors 2010 --------------------------------------------------
        % Driving function derived by approximation of the SDM
        %
        % 2.5D correction factor
        %        _______________________
        % g0 = \| -y_ref / (y_s - y_ref)
        %
        g0 = sqrt(- xref(1,2) / (xs(1,2) - xref(1,2)));
        %
        %                1 i w       ys    (2)/w       \
        % D_2.5D(x0,w) = - --- g0 ------- H1 | -|x0-xs| |
        %                2  c     |x0-xs|     \c       /
        %
        % See Spors and Ahrens (2010), eq.(24)
        %
        r = vector_norm(x0-xs,2);
        D = 1./2 .* 1i.*omega./c .* g0 ...
            .* xs(1,2)./r ...
            .* besselh(1,2,omega./c.*r);
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D point source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
