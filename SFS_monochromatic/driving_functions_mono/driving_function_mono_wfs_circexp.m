function D = driving_function_mono_wfs_circexp(x0,n0,Pm,mode,f,xq,conf)
%computes the wfs driving functions for a sound field expressed by cylindrical 
%expansion coefficients.
%
%   Usage: D = driving_function_mono_wfs_circexp(x0,n0,Pm,mode,f,xq,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       n0          - directions of the secondary sources / m [nx3]
%       Pm          - circular expansion coefficients of sound field
%       mode        - 'R' for regular expansion, 'S' for singular expansion
%       f           - frequency / Hz
%       xq          - optional expansion center coordinates, default: [0, 0, 0]
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_WFS_CIRCEXP(x0,n0,Pm,mode,f,xq,conf)
%
%   see also: driving_function_mono_wfs_sphexpS
%
%   References:
%     Spors2011 - "Local Sound Field Synthesis by Virtual Acoustic Scattering
%       and Time-Reversal" (AES131)
%     Gumerov,Duraiswami (2004) - "Fast Multipole Methods for the Helmholtz
%       Equation in three Dimensions", ELSEVIER

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
nargmin = 5;
nargmax = 7;
narginchk(nargmin,nargmax);
isargmatrix(x0,n0);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if nargin == nargmin
  xq = [0, 0, 0];
end
isargposition(xq);

%% ===== Configuration ==================================================
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;

%% ===== Computation ====================================================
Nce = (size(Pm, 1)-1)/2;

% apply shift to the center of expansion xq
x = x0(:,1)-xq(1);
y = x0(:,2)-xq(2);

% conversion to cylindrical coordinates
r0 = sqrt(x.^2 + y.^2);
phi0 = atan2(y,x);

% frequency depended stuff
omega = 2*pi*f;
k = omega/c;
kr0 = k.*r0;

% gradient in spherical coordinates
Gradr = zeros(size(x0,1),1);
Gradphi = zeros(size(x0,1),1);
Gradz = zeros(size(x0,1),1);

% directional weights for conversion spherical gradient into carthesian 
% coordinates + point product with normal vector n0 (directional derivative 
% in cartesian coordinates)
Sn0r     =  cos(phi0).*n0(:,1)...
         +  sin(phi0).*n0(:,2);
Sn0phi   = -sin(phi0).*n0(:,1)...
         +  cos(phi0).*n0(:,2);
Sn0z     =            n0(:,3);

% select suitable basis function
if strcmp('R', mode)
  circbasis = @besselj;
  circbasis_derived = @besselj_derived;
elseif strcmp('S', mode)
  circbasis = @(nu,z) besselh(nu,2,z);
  circbasis_derived = @(nu,z) besselh_derived(nu,2,z);
else
  error('unknown mode:');
end

%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% indexing the expansion coefficients
l = 0;

if strcmp('2D',dimension) || strcmp('3D',dimension)
  
  % === 2- or 3-Dimensional ============================================
  
  if (strcmp('default',driving_functions))    
    % --- SFS Toolbox ------------------------------------------------
    %                d    
    % D(x0, w) = -2 --- P(x0,w)
    %               d n
    % with cylindrical expansion of the sound field:
    %           \~~    oo       
    % P(x,w) =  >        B   F (x-xq) 
    %           /__ n=-oo  n  n
    %
    % where F = {R,S}.
    %
    % regular cylindrical basis functions:
    % 
    % R  (x) = J (kr)  . exp(j n phi))
    %  n        n       
    % singular cylindrical basis functions
    %          (2) 
    % S (x) = H   (kr) . exp(j n phi)
    %  n       n
 
    for n=-Nce:Nce
      l = l + 1;      
      cn_prime = k.*circbasis_derived(n,kr0);
      cn = circbasis(n,kr0);
      Yn = exp(1j.*n.*phi0);      
      Gradr   = Gradr   +       ( Pm(l).*cn_prime.*Yn  );
      Gradphi = Gradphi + 1./r0.*( Pm(l).*cn.*1j.*n.*Yn );
    end
    % directional gradient
    D = -2*( Sn0r.*Gradr + Sn0phi.*Gradphi + Sn0z.*Gradz );
  else
    error(['%s: %s, this type of driving function is not implemented ', ...
      'for 2D/3D.'],upper(mfilename),driving_functions);
  end  
else
  error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
