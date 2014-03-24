function D = driving_function_mono_nfchoa_ps(x0,xs,f,N,conf)
%DRIVING_FUNCTION_MONO_NFCHOA_PS returns the driving signal D for a point source
%in NFCHOA
%
%   Usage: D = driving_function_mono_nfchoa_ps(x0,xs,f,N,[conf])
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       xs          - position of virtual point source / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       N           - maximum order of spherical harmonics
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_PS(x0,xs,f,N,conf) returns NFCHOA driving
%   signals for the given secondary sources, the virtual point source position
%   and the frequency f.
%
%   References:
%       J. Ahrens (2012) - "Analytic Methods of Sound Field Synthesis", Springer
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, TU Berlin
%
%   see also: driving_function_mono_nfchoa, driving_function_imp_nfchoa_ps

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
isargmatrix(x0,xs);
isargpositivescalar(f,N);
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
X0 = conf.secondary_sources.center;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% secondary source positions
x00 = bsxfun(@minus,x0,X0);
[phi0,theta0,r0] = cart2sph(x00(:,1),x00(:,2),x00(:,3));
% point source
[phi,theta,r] = cart2sph(xs(:,1),xs(:,2),xs(:,3));
% wavenumber
omega = 2*pi*f;
% initialize empty driving signal
D = zeros(size(x0,1),1);

if strcmp('2D',dimension)
    
    % === 2-Dimensional ==================================================
    
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)
    
    % === 2.5-Dimensional ================================================
    
    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 2.5D point source, after Ahrens (2012), p.186 (5.8)
        %
        %                      __      (2)
        %               1     \       h|m| (w/c r)
        % D(phi0,w) = -----   /__    ------------- e^(i m (phi0-phi))
        %             2pi r0 m=-N..N  (2)
        %                             h|m| (w/c r0)
        %
        % see Wierstorf (2014), p.24 (2.41)
        for m=-N:N
            D = D + 1./(2.*pi.*r0) .* sphbesselh(abs(m),2,omega/c.*r) ./ ...
                sphbesselh(abs(m),2,omega/c.*r0) .* exp(1i.*m.*(phi0-phi));
        end
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)
    
    % === 3-Dimensional ==================================================
    
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % 3D point source, after Ahrens (2012), p.185 (5.7)
        %
        %                              N_   n_  (2)
        %                       1     \    \    hn (w/c r)   -m             
        % D(theta0,phi0,w) = -------  /__  /__  ------------ Yn(theta,phi) ...
        %                    2pi r0^2 n=0 m=-n  (2)
        %                                       hn (w/c r0)
        %                      m
        %                     Yn(theta0,phi0)
        %
        % see Wierstorf (2014), p.24 (2.40)
        %
        for n=0:N
            for m=-n:n
                D = D + 1./(2.*pi.*r0.^2) .* sphbesselh(n,2,omega/c.*r) ./ ...
                    sphbesselh(n,2,omega/c.*r0) .* ...
                    sphharmonics(n,-m,theta,phi) .* ...
                    sphharmonics(n,m,theta0,phi0);
            end
        end
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 3D point source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
