function D = driving_function_mono_nfchoa_pw(x0,nk,f,N,conf)
%DRIVING_FUNCTION_MONO_NFCHOA_PW driving signal for a plane wave in NFC-HOA
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
%   See also: driving_function_mono_nfchoa, driving_function_imp_nfchoa_pw

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
omega = 2*pi*f;

% modal window
win = modal_weighting(N,conf);

% Initialize empty driving signal
D = zeros(size(x0,1),1);

if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % 2D plane wave
        %
        %                      _N_
        %                2i    \       i^-m
        % D(phi0,w) = - -----  /__  ---------- e^(i m (phi0-phi_pw))
        %               pi r0  m=-N  (2)
        %                           Hm(w/c r0)
        %
        % https://sfs.rtfd.io/en/3.2/d_nfchoa/#equation-fd-nfchoa-plane-2d
        %
        for m=-N:N
            D = D - 2.*1i ./ (pi.*r0) ...
                .* win(abs(m)+1) ...
                .* (1i).^(-m) ./ besselh(m,2,omega./c.*r0) ...
                .* exp(1i.*m.*(phi0-phi_pw));
        end
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2D plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D plane wave
        %
        %                       _N_
        %                   2   \           i^-|m|
        % D_25D(phi0,w) = - --  /__  -------------------- e^(i m (phi0-phi_pw))
        %                   r0  m=-N       (2)
        %                             i w/c h|m|(w/c r0)
        %
        % https://sfs.rtfd.io/en/3.2/d_nfchoa/#equation-fd-nfchoa-plane-25d
        %
        for m=-N:N
            D = D - 2./r0 ...
                .* win(abs(m)+1) ...
                .* (1i).^(-abs(m)) ...
                ./ (1i .* omega/c .* sphbesselh(abs(m),2,omega./c.*r0)) ...
                .* exp(1i.*m.*(phi0-phi_pw));
        end
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % 3D plane wave
        %
        %                            _N_  _n_         -m
        %                      4pi   \    \   i^(-n) Yn(theta_pw,phi_pw)
        % D(theta0,phi0,w) = - ----  /__  /__ -------------------------- ...
        %                      r0^2  n=0 m=-n            (2)
        %                                         i w/c hn(w/c r0)
        %                       m
        %                    x Yn(theta0,phi0)
        %
        % https://sfs.rtfd.io/en/3.2/d_nfchoa/#equation-fd-nfchoa-plane-3d
        %
        for n=0:N
            for m=-n:n
                D = D - 4.*pi ./ r0.^2 ...
                    .* win(n+1) ...
                    .* (1i).^(-n).*sphharmonics(n,-m,theta_pw,phi_pw) ...
                    ./ (1i .* omega./c .* sphbesselh(n,2,omega./c.*r0)) ...
                    .* sphharmonics(n,m,theta0,phi0);
            end
        end
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 3D plane wave.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
