function G = greens_function_mono(x,y,z,xs,src,f,conf)
%GREENS_FUNCTION_MONO returns a Green's function in the frequency domain
%
%   Usage: G = greens_function_mono(x,y,z,xs,src,f,conf)
%
%   Input options:
%       x,y,z   - x,y,z points for which the Green's function should be
%                 calculated / m
%       xs      - position of the source
%       src     - source model of the Green's function. Valid models are:
%                   'ps'  - point source
%                   'ls'  - line source
%                   'pw'  - plane wave
%                   'dps' - dipole point source
%       f       - frequency of the source / Hz
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       G       - Green's function evaluated at the points x,y,z
%
%   GREENS_FUNCTION_MONO(x,y,z,xs,src,f,conf) calculates the Green's function
%   for the given source model located at xs for the given points x,y and the
%   frequency f.
%
%   References:
%       H. Wierstorf, J. Ahrens, F. Winter, F. Schultz, S. Spors (2015) -
%       "Theory of Sound Field Synthesis"
%       J. Ahrens (2012) - "Analytic Methods of Sound Field Synthesis"
%
%   See also: sound_field_mono

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
% Disabled checking for performance reasons


%% ===== Configuration ==================================================
c = conf.c;
phase = conf.phase;


%% ===== Computation =====================================================
% Frequency
omega = 2*pi*f;
% Calculate Green's function for the given source model
if strcmp('ps',src)
    % Source model for a point source: 3D Green's function.
    %
    %              1  e^(-i w/c |x-xs|)
    % G(x-xs,w) = --- -----------------
    %             4pi      |x-xs|
    %
    % see: Wierstorf et al. (2015), eq.(#S:ps)
    %
    G = 1/(4*pi) * exp(-1i*omega/c .* sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2)) ./ ...
            sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2);

elseif strcmp('dps',src)
    % Source model for a dipole point source: derivative of 3D Green's function.
    %
    %  d                1   / iw       1    \   (x-xs) ns
    % ---- G(x-xs,w) = --- | ----- + ------- | ----------- e^(-i w/c |x-xs|)
    % d ns             4pi  \  c     |x-xs| /   |x-xs|^2
    %
    % see Wierstorf et al. (2015), eq.(#S:dps)
    %
    % r = |x-xs|
    r = sqrt((x-xs(1)).^2+(y-xs(2)).^2+(z-xs(3)).^2);
    % scalar = (x-xs) nxs
    scalar = xs(4).*(x-xs(1)) + xs(5).*(y-xs(2))  + xs(6).*(z-xs(3));
    %
    G = 1/(4*pi) .* (1i*omega/c + 1./r) .* scalar./r.^2 .* exp(-1i*omega/c.*r);

elseif strcmp('mps',src)
    % Source model for a uniformly moving point source: retarded 3D Green's 
    % function.
    %
    %                 1   e^(+-i w/c R)
    % G(x-xs,vs,w) = --- --------------
    %                4pi      R_1
    %
    % vs = v * ns   % velocity vector
    % M = v/c  % Mach number
    % xparallel = (x-xs)*ns  % component in direction of movement
    % R_1 = sqrt( M^2*xorth.^2 + (1-M^2)*|x-xs|.^2 )
    % R = (M*xorth + R1)./(1-M*2);
    %    
    % see: Ahrens (2012), eq.(5.60)
    %
    v = norm(xs(4:6));  % velocity of sound source
    nxs = xs(4:6)./v;  % direction of movement
    xs = xs(1:3);  % 
    M = v/c;  % 
    % shift coordinates
    x = x-xs(1);
    y = y-xs(2);
    z = z-xs(3);
    % component of x in direciton  of movement: scalar = x*nxs  
    xparallel = nxs(1).*x + nxs(2).*y  + nxs(3).*z;
    
    R1 = sqrt( M^2.*xparallel.^2 + (1-M^2)*(x.^2 + y.^2 + z.^2) );
    Rplus = (M*xparallel + R1)./(1-M.^2);  
    G = 1/(4*pi) * exp(-1i*omega/c.*Rplus)./R1;
elseif strcmp('ls',src)
    % Source model for a line source: 2D Green's function.
    %
    %                i   (2) / w        \
    % G(x-xs,w) =  - -  H0  |  - |x-xs|  |
    %                4       \ c        /
    %
    % see: Wierstorf et al. (2015), eq.(#S:ls)
    %
    G = -1i/4 * besselh(0,2,omega/c* ...
        sqrt( (x-xs(1)).^2 + (y-xs(2)).^2 + (z-xs(3)).^2 ));

elseif strcmp('pw',src)
    % Source model for a plane wave:
    %
    % G(x,w) = e^(-i w/c n x)
    %
    % see: Wierstorf et al. (2015), eq.(#S:pw)
    %
    % Direction of plane wave
    nxs = xs(:,1:3) / norm(xs(:,1:3));
    %
    % The following code enables us to replace this two for-loops
    % for ii = 1:size(x,1)
    %     for jj = 1:size(x,2)
    %         S(ii,jj) = exp(-1i*omega/c.*nxs*[x(ii,jj) y(ii,jj) z(ii,jj)]');
    %     end
    % end
    %
    % Get a matrix in the form of
    % 1 1 1 0 0 0 0 0 0
    % 0 0 0 1 1 1 0 0 0
    % 0 0 0 0 0 0 1 1 1
    E = eye(3*size(x,1));
    E = E(1:3:end,:)+E(2:3:end,:)+E(3:3:end,:);
    % Multiply this matrix with the plane wave direction
    N = repmat(nxs,size(x,1)) .* E;
    % Interlace x,y,z into one matrix
    % x11 x12 ... x1m
    % y11 y12 ... y1m
    % z11 z12 ... z1m
    % .   .       .
    % .   .       .
    % xn1 xn2 ... xnm
    % yn1 yn2 ... ynm
    % zn1 zn2 ... znm
    XYZ = zeros(3*size(x,1),size(x,2));
    XYZ(1:3:end,:) = x;
    XYZ(2:3:end,:) = y;
    XYZ(3:3:end,:) = z;
    % Calculate sound field
    G = exp(-1i*omega/c.*N*XYZ);

else
    error('%s: %s is not a valid source model for the Green''s function', ...
        upper(mfilename),src);
end

% Add phase to be able to simulate different time steps
G = G .* exp(-1i*phase);
