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
%       H. Wierstorf, J. Ahrens, F. Winter, F. Schultz, S. Spors (2015) -
%       "Theory of Sound Field Synthesis"
%       S. Spors, R. Rabenstein, J. Ahrens (2008) - "The Theory of Wave Field
%       Synthesis Revisited", AES124
%       E. Start (1997) - "Direct Sound Enhancement by Wave Field Synthesis", 
%       PhD thesis, TU Delft
%       E. Verheijen (1997) - "Sound Reproduction by Wave Field Synthesis", PhD
%       thesis, TU Delft
%       S. Spors, J. Ahrens (2010) - "Analysis and Improvement of
%       Pre-equalization in 2.5-Dimensional Wave Field Synthesis", AES128
%
%   See also: driving_function_mono_wfs, driving_function_imp_wfs_ps

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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

% Frequency
omega = 2*pi*f;


if strcmp('2D',dimension) || strcmp('3D',dimension)

    % === 2- or 3-Dimensional ============================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % D using a point sink and large distance approximation
        %
        %            1  i w  (x0-xs) nx0
        % D(x0,w) = --- --- ------------- e^(-i w/c |x0-xs|)
        %           2pi  c  |x0-xs|^(3/2)
        %
        % see Wierstorf et al. (2015), eq.(#D:wfs:ps)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = 1/(2*pi) .* (1i*omega)/c .* ...
            vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(-1i*omega/c.*r);
        %
    elseif strcmp('point_source',driving_functions)
        % D using a point source as source model
        %
        %            1  / i w      1    \  (x0-xs) nx0
        % D(x0,w) = --- | --- - ------- |  ----------- e^(-i w/c |x0-xs|)
        %           2pi \  c    |x0-xs| /   |x0-xs|^2
        %
        % see Wierstorf et al. (2015), eq.(#D:wfs:ps:woapprox)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = 1/(2*pi) .* ( (1i*omega)/c - 1./r ) .* ...
            vector_product(x0-xs,nx0,2) ./ r.^2 .* exp(-1i*omega/c.*r);
        %
    elseif strcmp('line_source',driving_functions)
        % D using a line source as source model (truly 2D model)
        %
        %             1  i w (x0-xs) nx0  (2)/ w         \
        % D(x0,x) = - -- --- ----------- H1  | - |x0-xs| |
        %             2c  c    |x0-xs|       \ c         /
        %
        % see Spors et al. (2008), eq.(23)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Driving signal
        D = -1/(2*c) .* 1i*omega/c * vector_product(x0-xs,nx0,2) ./ r .* besselh(1,2,omega/c*r);
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
      % see Start (1997), eq. (3.11)
      %
      g0 = sqrt( vector_norm(xref-x0,2) ./ (vector_norm(xref-x0,2) + r) );
      %
      % D_2.5D(x0,w) =
      %       ___    ___
      %      | 1    |i w (x0-xs) nx0
      % g0 _ |--- _ |--- ------------- e^(-i w/c |x0-xs|)
      %     \|2pi  \| c  |x0-xs|^(3/2)
      %
      % see Start (1997), eq. (3.10)
      %
      % Driving signal
      D = g0/sqrt(2*pi) .* sqrt( 1i*omega/c ) .* ...
        vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(-1i*omega/c.*r);
      %
    case {'reference_line', 'delft1988'}
      % Driving function with two stationary phase approximations,
      % reference to a line parallel to a LINEAR secondary source distribution
      %
      % distance ref-line to linear ssd
      dref = vector_product(xref-x0,nx0,2);
      % distance source and linear ssd
      ds = vector_product(xs-x0,nx0,2);
      %
      % 2.5D correction factor
      %        _______________________
      % g0 = \| d_ref / (d_ref - d_s)
      %
      % see Start (1997), eq. (3.16)
      %
      g0 = sqrt( dref ./ (dref - ds) );
      %                      _____
      %                     |i w   (x0-xs) nx0
      % D_2.5D(x0,w) = g0 _ |----- ------------  e^(-i w/c |x0-xs|)
      %                    \|2pi c |x0-xs|^(3/2)
      %
      % see Start (1997), eq. (3.17)
      %
      % r = |x0-xs|
      r = vector_norm(x0-xs,2);
      % Driving signal
      D = sqrt(1i*omega/(2*pi*c)) .* g0 .* ...
        vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(-1i*omega/c .* r);
      %
    case 'legacy'
      % --- SFS Toolbox ------------------------------------------------
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
      % see Wierstorf et al. (2015), eq.(#D:wfs:ps:2.5D)
      %
      % r = |x0-xs|
      r = vector_norm(x0-xs,2);
      % Driving signal
      D = g0/(2*pi) .* sqrt(1i*omega/c) .* ...
        vector_product(x0-xs,nx0,2) ./ r.^(3/2) .* exp(-1i*omega/c.*r);
      %
    case 'point_source'
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
      % see Wierstorf et al. (2015), eq.(#D:wfs:ps:woapprox:2.5D)
      %
      % r = |x0-xs|
      r = vector_norm(x0-xs,2);
      % Driving signal
      D = g0/(2*pi) .* ( sqrt(1i*omega/c) - sqrt(c/(1i*omega) ./ r ) ) .* ...
        vector_product(x0-xs,nx0,2) ./ r.^2 .* exp(-1i*omega/c .* r);
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
      % see Spors and Ahrens (2010), eq.(24)
      %
      r = vector_norm(x0-xs,2);
      D = 1/2 * 1i*omega/c * g0 * xs(1,2)./r .* besselh(1,2,omega/c*r);
      %
    otherwise
      error(['%s: %s, this type of driving function is not implemented ', ...
        'for a 2.5D point source.'],upper(mfilename),driving_functions);
  end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
