function D = driving_function_mono_nfchoa_circexp(x0, Pm,f,conf)
%computes the nfchoa driving functions for a sound field expressed by regular
%spherical expansion coefficients.
%
%   Usage: D = driving_function_mono_nfchoa_circexp(x0,Pm,f,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [N0x3]
%       Pm          - regular circular expansion coefficients of sound field
%                     [N x Nf]
%       f           - frequency / Hz [Nf x 1] or [1 x Nf]
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_CIRCEXP(x0,Pm,f,conf) returns the NFCHOA
%   driving signals for the given secondary sources x0, the virtual sound ex-
%   pressed by regular circular expansion coefficients Pm, and the frequency f.
%
%   see also: driving_function_mono_wfs_circexp
%             driving_function_mono_nfchoa_sphexp

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(Pm,x0);
isargvector(f);
isargstruct(conf);

%% ===== Configuration ==================================================
c = conf.c;
dimension = conf.dimension;
Xc = conf.secondary_sources.center;

%% ===== Computation ====================================================
Nce = (size(Pm, 1)-1)/2;
% Calculate the driving function in time-frequency domain

% secondary source positions
x00 = bsxfun(@minus,x0(:,1:3),Xc);
[phi0, ~,r0] = cart2sph(x00(:,1),x00(:,2),x00(:,3));

% frequency depended stuff
omega = 2*pi*row_vector(f);  % [1 x Nf]
k = omega./c;  % [1 x Nf]
kr0 = r0 * k;  % [N0 x Nf]

% initialize empty driving signal
D = zeros(size(kr0));  % [N0 x Nf]

% indexing the expansion coefficients
l = 0;

% Calculate the driving function in time-frequency domain

switch dimension
  case '2D'
    % === 2-Dimensional ==================================================
    %
    %                       __      P
    %               1      \         m
    % D(phi0,w) = ------   /__    -------- e^(i m phi0)
    %             2pi r0  m=-N..N   G
    %                                m
    %
    % with regular circular expansion of the sound field:
    %          \~~ oo                 
    % P(x,w) =  >        P  J (kr) . e^(+j m phi)
    %          /__ m=-oo  m  m
    %
    % and 3D free field Green's Function:
    %            \~~ oo                 
    % G  (x,w) =  >        G  J (kr) . e^(+j m phi)
    %  ls        /__ m=-oo  m  m
    %
    % with the regular expansion coefficients:
    %    m     -i      (2)
    %   G  =  ----- . H   (kr0)
    %    n      4      m
    for m=-Nce:Nce
      l = l+1;
      D = D + bsxfun(@times, Pm(l,:), exp(1i.*m.*phi0))./ besselh(m,2,kr0);
    end
    % factor from expansion of 2D free field Green's Function
    D = D./(-1j/4);
    % normalization due to size of circular array
    D = bsxfun(@rdivide, D, 2*pi*r0);
  case {'2.5D', '3D'}
    % convert circular expansion to spherical expansion
    Pmn = sphexp_convert_circexp(Pm);
    % compute driving function for spherical expansion
    D = driving_function_mono_nfchoa_sphexp(x0,Pmn,f,conf);
  otherwise
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
