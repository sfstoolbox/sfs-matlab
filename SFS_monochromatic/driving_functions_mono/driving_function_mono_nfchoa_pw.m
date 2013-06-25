function D = driving_function_mono_nfchoa_pw(x0,nk,f,conf)
%DRIVING_FUNCTION_MONO_NFCHOA_PW returns the driving signal D for a plane wave
%in NFCHOA
%
%   Usage: D = driving_function_mono_nfchoa_pw(x0,nk,f,[conf])
%
%   Input parameters:
%       x0          - position of the secondary sources (m) [nx3]
%       nk          - direction of virtual plane wave (m) [nx3]
%       f           - frequency of the monochromatic source (Hz)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_NFCHOA_PW(x0,nk,f,src,conf) returns NFCHOA driving
%   signals for the given secondary sources, the virtual plane wave direction
%   and the frequency f.
%
%   References:
%       FIXME: Update references
%
%   see also: driving_function_mono_nfchoa, driving_function_imp_nfchoa_pw

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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(x0,nk);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
xref = conf.xref;
X0 = conf.X0;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% angle of the secondary sources
[alpha_x0,beta_x0] = cart2sph(bsxfun(@minus,x0,X0));
% angle of plane wave
[alpha_pw,beta_pw] = cart2sph(nk);
% wavenumber
k = 2*pi*f/c;
% max order of spherical harmonics
M = size(x0,1);
if(isodd(M))
    N=(M-1)/2;
else
    N=floor((M+1)/2);
end
% initialize empty driving signal
D = zeros(M,1);

if strcmp('2D',dimension)
    
    % === 2-Dimensional ==================================================
    
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2D plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)
    
    % === 2.5-Dimensional ================================================
    
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        %
        %       __           4pi i^|n|
        %      \      --------------------- e^(i n (alpha_x0-alpha_pw) )
        % D =  /__          (2) 
        %    n=-N..N  -ik H|n| (k|x0-xref|)  
        %                      
        % R = |x0-xref|
        % NOTE: it makes only sense to use the center point as reference point.
        % Otherwise we will have no radius at all.
        R = norm(x0(1,:)-X0,2);
        for n=-N:N
            D = D + 4.*pi .* 1i.^(-abs(n)) ./ ...
                ( -1i .* k .* sphbesselh(abs(n),2,k.*R) ) .* exp(1i.*n.*(alpha_x0-alpha_pw));
        end
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)
    
    % === 3-Dimensional ==================================================
    
    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 3D plane wave.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
