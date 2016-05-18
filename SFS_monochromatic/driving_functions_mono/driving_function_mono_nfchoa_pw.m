function D = driving_function_mono_nfchoa_pw(x0,nk,f,N,conf)
%DRIVING_FUNCTION_MONO_NFCHOA_PW returns the driving signal D for a plane wave
%in NFCHOA
%
%   Usage: D = driving_function_mono_nfchoa_pw(x0,nk,f,N,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nk          - direction of virtual plane wave / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       N           - maximum order of spherical harmonics
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_PW(x0,nk,f,N,conf) returns NFCHOA driving
%   signals for the given secondary sources, the virtual plane wave direction
%   and the frequency f.
%
%   See also: driving_function_mono_nfchoa, driving_function_imp_nfchoa_pw

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
isargmatrix(x0,nk);
isargpositivescalar(f,N);
isargstruct(conf);


%% ===== Configuration ==================================================
X0 = conf.secondary_sources.center;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation ====================================================
% Angle of the secondary sources
points = bsxfun(@minus,x0,X0);
[phi0,theta0,r0] = cart2sph(points(:,1),points(:,2),points(:,3));

% Angle of plane wave
[phi_pw,theta_pw,~] = cart2sph(nk(:,1),nk(:,2),nk(:,3));

% Wavenumber
w = 2*pi*f;

% Initialize empty driving signal
D = zeros(size(x0,1),1);

if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 2D plane wave, see
        % http://sfstoolbox.org/doc/latest/math/#equation-D.nfchoa.pw.2D
        %
        %                      _N_
        %                2i    \       i^-m
        % D(phi0,w) = - -----  /__  ---------- e^(i m (phi0-phi_pw))
        %               pi r0  m=-N  (2)
        %                           Hm(w/c r0)
        %
        for m=-N:N
            D = D - 2.*1i./(pi.*r0) .* 1i^(-m)./besselh(m,2,w/c.*r0) .* ...
                exp(1i.*m.*(phi0-phi_pw));
        end
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2D plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D plane wave, see
        % http://sfstoolbox.org/doc/latest/math/#equation-D.nfchoa.pw.2.5D
        %
        %                     _N_
        %                 2i  \          i^-|m|
        % D_25D(phi0,w) = --  /__  ------------------ e^(i m (phi0-phi_pw) )
        %                 r0  m=-N       (2)
        %                           w/c h|m|(w/c r0)
        %
        for m=-N:N
            D = D + 2.*1i./r0 .* 1i.^(-abs(m)) ./ ...
                ( w/c .* sphbesselh(abs(m),2,w/c.*r0) ) .* ...
                exp(1i.*m.*(phi0-phi_pw));
        end
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 3D plane wave, see
        % http://sfstoolbox.org/doc/latest/math/#equation-D.nfchoa.pw.3D
        %
        %                          _N_  _n_         -m
        %                     4i   \    \   i^(-n) Yn(theta_pw,phi_pw)
        % D(theta0,phi0,w) = ----  /__  /__ -------------------------- ...
        %                    r0^2  n=0 m=-n           (2)
        %                                        w/c hn(w/c r0)
        %                       m
        %                    x Yn(theta0,phi0)
        %
        for n=0:N
            for m=-n:n
                D = D + 4.*1i./r0.^2 .* 1i.^(-n).*sphharmonics(n,-m,theta_pw,phi_pw) ./...
                    ( w./c .* sphbesselh(n,2,w./c.*r0) ) .* ...
                    sphharmonics(n,m,theta0,phi0);
            end
        end
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 3D plane wave.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
