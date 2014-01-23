function D = driving_function_mono_wfs_ps(x0,nx0,xs,f,conf)
%DRIVING_FUNCTION_MONO_WFS_PS returns the driving signal D for a point source in
%WFS
%
%   Usage: D = driving_function_mono_wfs_ps(x0,nx0,xs,f,[conf])
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nx0         - directions of the secondary sources / m [nx3]
%       xs          - position of virtual point source / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_WFS_PS(x0,xs,f,src,conf) returns WFS driving signals
%   for the given secondary sources, the virtual point source position and the
%   frequency f.
%
%   References:
%       FIXME: Update references
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%           2.5-Dimensional Wave Field Synthesis (AES128)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: driving_function_mono_wfs, driving_function_imp_wfs_ps

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
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
        % D using a point sink and large distance approximation
        %
        %              1  i w  (x0-xs) nx0
        % D(x0,w) = - --- --- ------------- e^(-i w/c |x0-xs|)
        %             2pi  c  |x0-xs|^(3/2)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -1/(2*pi) .* (1i*omega)/c .* ...
            vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(-1i*omega/c.*r);
        %
    elseif strcmp('point_source',driving_functions)
        % D using a point source as source model
        %
        %              1  / i w      1    \  (x0-xs) nx0
        % D(x0,w) = - --- | --- - ------- |  ----------- e^(-i w/c |x0-xs|)
        %             2pi \  c    |x0-xs| /   |x0-xs|^2
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -1/(2*pi) .* ( (1i*omega)/c - 1./r ) .* ...
            vector_product(x0-xs,nx0,2) ./ r.^2 .* exp(-1i*omega/c.*r);
        %
    elseif strcmp('line_source',driving_functions)
        % D using a line source as source model (truly 2D model)
        % Spors et al, The theory of Wave Field Synthesis Revisited, 2008,
        % AES. Equation (23)
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = 1i*omega/c * vector_product(x0-xs,nx0,2) ./ r .* besselh(1,2,omega/c*r);
        %
        %
    elseif strcmp('delft1988',driving_functions)
        % --- Delft 1988 -------------------------------------------------
        to_be_implemented;
        %
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)
    
    % === 2.5-Dimensional ================================================
    
    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D correction factor
        %        _____________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using a point source and large distance approximation
        %                         ___
        %                  g0    |i w  (x0-xs) nx0
        % D_2.5D(x0,w) = - --- _ |--- ------------- e^(-i w/c |x0-xs|)
        %                  2pi  \| c  |x0-xs|^(3/2)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = -g0/(2*pi) .* sqrt(1i*omega/c) .* ...
            vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(-1i*omega/c.*r);
        %
    elseif strcmp('point_source',driving_functions)
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % D_2.5D using a point source as source model
        %
        % D_2.5D(x0,w) =
        %             ___       ___
        %    g0  /   |i w      | c      1    \  (x0-xs) nx0
        %   ---  | _ |---  - _ |---  ------- |  ----------- e^(-i w/c |x0-xs|)
        %   2pi  \  \| c      \|i w  |x0-xs| /   |x0-xs|^2
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = g0/(2*pi) .* ( sqrt(1i*omega/c) - sqrt(c/(1i*omega) ./ r ) ) .* ...
            vector_product(x0-xs,nx0,2) ./ r.^2 .* exp(-1i*omega/c .* r);
        %
    elseif strcmp('delft1988',driving_functions)
        % --- Delft 1988 -------------------------------------------------
        % Verheijen, Sound Reproduction by Wave Field Synthesis, PhD Thesis
        % Equation (2.27)
        % Note: So far this works only for a linear secondary source distribution
        %       located on the y-axis
        %
        % 2.5D correction factor
        %        _______________________
        % g0 = \| -y_ref / (y_s - y_ref)
        %
        g0 = sqrt(- xref(1,2) / (xs(1,2) - xref(1,2)));
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = sqrt(1i*omega/c/(2*pi)) * g0 * vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(-1i*omega/c .* r);
        %
    elseif strcmp('opperschall',driving_functions)
        % --- Opperschall -------------------------------------------------
        % Opperschall, Equation (3.14)
        % Note: Driving function with only one stationary phase
        % approximation, reference to one point in field
        %
        % 2.5D correction factor
        g0 = sqrt( vector_norm(x0-xref,2) ./ (vector_norm(xs-x0,2) + vector_norm(x0-xref,2)) );
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % driving signal
        D = sqrt(2*pi*1i*omega/c) * g0 .* vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(-1i*omega/c .* r);
        %
    elseif strcmp('volk2010',driving_functions)
        % --- Voelk 2010 --------------------------------------------------
        %
        % D_2.5D(x0,w) = 
        %    ___    ___    _____________________
        %   | 1    |i w   |      |xref-x0|         (x0-xs) nx0
        % _ |--- _ |--- _ |---------------------  ------------- e^(-i w/c |x0-xs|)
        %  \|2pi  \| c   \| |x0-xs| + |xref-x0|   |x0-xs|^(3/2)
        %
        g0 = vector_norm(xref-x0,2);
        r = vector_norm(x0-xs,2);
        D = 1/sqrt(2*pi) * sqrt(1i*omega/c) * sqrt(g0/(r+g0)) * ...
            vector_product(x0-xs,nx0,2)./r.^(3/2) .* exp(-1i*omega/c.*r); 
        %
    elseif strcmp('SDMapprox',driving_functions)
        % --- Spors 2010 --------------------------------------------------
        % Driving function derived by approximation of the SDM, eq(24) in
        % AES Paper "Analysis and Improvement of Pre-Equalization in 2.5-
        % Dimensional Wave Field Synthesis
        g0 = sqrt(- xref(1,2) / (xs(1,2) - xref(1,2)));
        r = vector_norm(x0-xs,2);
        D = 1i*omega/c * g0 * xs(1,2)./r .* besselh(1,2,omega/c*r);
        %
        %
    
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D point source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
