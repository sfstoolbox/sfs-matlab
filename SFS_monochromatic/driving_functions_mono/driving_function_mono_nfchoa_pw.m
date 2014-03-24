function D = driving_function_mono_nfchoa_pw(x0,nk,f,N,conf)
%DRIVING_FUNCTION_MONO_NFCHOA_PW returns the driving signal D for a plane wave
%in NFCHOA
%
%   Usage: D = driving_function_mono_nfchoa_pw(x0,nk,f,N,[conf])
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nk          - direction of virtual plane wave / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       N           - maximum order of spherical harmonics
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_PW(x0,nk,f,N,conf) returns NFCHOA driving
%   signals for the given secondary sources, the virtual plane wave direction
%   and the frequency f.
%
%   References:
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, TU Berlin
%
%   see also: driving_function_mono_nfchoa, driving_function_imp_nfchoa_pw

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
isargmatrix(x0,nk);
isargpositivescalar(f,N);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
X0 = conf.secondary_sources.center;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% angle of the secondary sources
points = bsxfun(@minus,x0,X0);
[phi0,theta0,r0] = cart2sph(points(:,1),points(:,2),points(:,3));
% angle of plane wave
[phi_pw,theta_pw,~] = cart2sph(nk(:,1),nk(:,2),nk(:,3));
% wavenumber
w = 2*pi*f;
% initialize empty driving signal
D = zeros(size(x0,1),1);

if strcmp('2D',dimension)
    
    % === 2-Dimensional ==================================================
    
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        %                        __
        %                2i     \        i^-m
        % D(phi0,w) = - -----   /__   ----------  e^(i m (phi0-phi_pw))
        %               pi r0 m=-N..N  (2)
        %                             Hm  (w/c r0)
        %
        % see Wierstorf (2014), p.23 (2.36)
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
        %                      __
        %                 2i  \            i^|m|
        % D_25D(phi0,w) = --  /__    ------------------ e^(i m (phi0-phi_pw) )
        %                 r0 m=-N..N       (2)
        %                             w/c h|m| (w/c r0)  
        % 
        % see wierstorf (2014), p.23 (2.37)
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
        %                           __    __             -m
        %                    2i    \     \       i^(-n) Yn (theta_pw,phi_pw)  m
        % D(theta0,phi0,w) = ----  /__   /__     --------------------------- Yn (theta0,phi0)
        %                    r0^2 n=0..N m=-n..n           (2)
        %                                             w/c hn  (w/c r0)
        %
        % see Wierstorf (2014), p.23 (2.35)
        %
        for n=0:N
            for m=-n:n
                D = D + 2.*1i./r0.^2 .* 1i.^(-n).*sphharmonics(n,-m,theta_pw,phi_pw) ./...
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
