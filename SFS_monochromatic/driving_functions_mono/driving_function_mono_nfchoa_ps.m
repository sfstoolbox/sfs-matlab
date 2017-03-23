function D = driving_function_mono_nfchoa_ps(x0,xs,f,N,conf)
%DRIVING_FUNCTION_MONO_NFCHOA_PS returns the driving signal D for a point source
%in NFCHOA
%
%   Usage: D = driving_function_mono_nfchoa_ps(x0,xs,f,N,conf)
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       xs          - position of virtual point source / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       N           - maximum order of spherical harmonics
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_PS(x0,xs,f,N,conf) returns NFCHOA driving
%   signals for the given secondary sources, the virtual point source position
%   and the frequency f.
%
%   See also: driving_function_mono_nfchoa, driving_function_imp_nfchoa_ps

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
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
isargmatrix(x0,xs);
isargpositivescalar(f,N);
isargstruct(conf);


%% ===== Configuration ==================================================
xref = conf.xref;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;
X0 = conf.secondary_sources.center;


%% ===== Computation ====================================================

% Secondary source positions
x00 = bsxfun(@minus,x0,X0);
[phi0,theta0,r0] = cart2sph(x00(:,1),x00(:,2),x00(:,3));

% Point source
xs0 = bsxfun(@minus,xs,X0);
[phi,theta,r] = cart2sph(xs0(:,1),xs0(:,2),xs0(:,3));

% Wavenumber
omega = 2*pi*f;

% Initialize empty driving signal
D = zeros(size(x0,1),1);

if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D point source
        %
        %                     _N_    (2)
        %               1     \     h|m|(w/c r)
        % D(phi0,w) = ------  /__  ------------- e^(i m (phi0-phi))
        %             2pi r0  m=-N   (2)
        %                           h|m|(w/c r0)
        %
        % See http://sfstoolbox.org/#equation-D.nfchoa.ps.2.5D
        for m=-N:N
            D = D + 1 ./ (2.*pi.*r0) ...
                .* sphbesselh(abs(m),2,omega./c.*r) ...
                ./ sphbesselh(abs(m),2,omega./c.*r0) ...
                .* exp(1i.*m.*(phi0-phi));
        end
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % 3D point source
        %
        %                              _N_  _n_   (2)
        %                       1      \    \    hn(w/c r)   -m
        % D(theta0,phi0,w) = -------   /__  /__ ----------- Yn(theta,phi) ...
        %                    2pi r0^2  n=0 m=-n   (2)
        %                                        hn(w/c r0)
        %                       m
        %                    x Yn(theta0,phi0)
        %
        % See http://sfstoolbox.org/#equation-D.nfchoa.ps.3D
        %
        for n=0:N
            for m=-n:n
                D = D + 1 ./ (2.*pi.*r0.^2) ...
                    .* sphbesselh(n,2,omega./c.*r) ...
                    ./ sphbesselh(n,2,omega./c.*r0) ...
                    .* sphharmonics(n,-m,theta,phi) ...
                    .* sphharmonics(n,m,theta0,phi0);
            end
        end
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 3D point source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
