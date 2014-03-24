function D = driving_function_mono_wfs_fs(x0,nx0,xs,f,conf)
%DRIVING_FUNCTION_MONO_WFS_FS returns the driving signal D for a focused source
%in WFS
%
%   Usage: D = driving_function_mono_wfs_fs(x0,nx0,xs,f,[conf])
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nx0         - directions of the secondary sources / m [nx3]
%       xs          - position of virtual focused source / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_WFS_FS(x0,xs,f,src,conf) returns WFS driving signals
%   for the given secondary sources, the virtual focused source position and the
%   frequency f.
%
%   References:
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, TU Berlin
%       S. Spors, H. Wierstorf, M. Geier, J. Ahrens (2009) - "Physical and
%       Perceptual Properties of Focused Sources in Wave Field Synthesis",
%       AES127
%
%   see also: driving_function_mono_wfs, driving_function_imp_wfs_ps

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargmatrix(x0,nx0,xs);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
xref = conf.xref;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% frequency
omega = 2*pi*f;


if strcmp('2D',dimension) || strcmp('3D',dimension)
    
    % === 2- or 3-Dimensional ============================================
    
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % D using an approximated point sink as source model
        %                  
        %              1  i w (x0-xs) nx0
        % D(x0,w) = - --- --- ------------- e^(i w/c |x0-xs|)
        %             2pi  c  |x0-xs|^(3/2)
        %
        % see Wierstorf (2014), p. 27 (2.62)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -1/(2*pi) * 1i*omega/c .* ...
            vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(1i*omega/c.*r);
        %
    elseif strcmp('point_sink',driving_functions)
        % D using a point sink as source model
        %
        % D(x0,w) =
        %                    
        %    1  /i w      1    \  (x0-xs) nx0
        % - --- |--- + ------- |  ----------- e^(i w/c |x0-xs|)
        %   2pi \ c    |x0-xs| /   |x0-xs|^2
        %
        % see Wierstorf (2014), p. 27 (2.60)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -1/(2*pi) * ( -1i.*omega/c + 1./r ) .* ...
            vector_product(x0-xs,nx0,2) ./ r.^2 .* exp(1i*omega./c.*r);
        %
    elseif strcmp('line_sink',driving_functions)
        % D using a line sink as source model
        %
        %              iw (x0-xs)nk  (1)/ w         \
        % D(x0,w) =  - -- --------- H1  | - |x0-xs| |
        %              2c  |x0-xs|      \ c         /
        %
        % see Wierstorf (2014), p. 27 (2.66)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
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
        %                -g0    |i w (x0-xs) nx0
        % D_2.5D(x0,w) = ---  _ |--- ------------- e^(i w/c |x0-xs|)
        %                2pi   \| c  |x0-xs|^(3/2)
        %
        % see Wierstorf (2014), p.27 (2.63)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -g0/(2*pi) * sqrt( 1i*omega/c ) .* ...
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
        % see Spors (2009), (7)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
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
        % see Spors et al. (2009), (6)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
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
        %   -g0  /   |i w|     | c |    1    \  (x0-xs) nx0
        %   ---  | _ |---  + _ |---  ------- |  ----------- e^(i w/c |x0-xs|)
        %   2pi  \  \| c      \|i w  |x0-xs| /   |x0-xs|^2
        %
        % see Wierstorf (2014), p.27 (2.61)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -g0/(2*pi) .* ( sqrt(1i*omega/c) + sqrt(c/(1i*omega)) ./ r ) .* ...
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
