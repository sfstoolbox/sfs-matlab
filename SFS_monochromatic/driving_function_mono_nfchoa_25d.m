function [D] = driving_function_mono_nfchoa_25d(x0,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_NFCHOA_25D returns the driving signal D for 2.5D HOA
%
%   Usage: D = driving_function_mono_nfchoa_25d(x0,xs,src,f,[conf])
%
%   Input parameters:
%       x0          - position and direction of the secondary sources (m)
%       xs          - position of virtual source or direction of plane wave (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       f           - frequency of the monochromatic source (Hz)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal (1x1)
%
%   DRIVING_FUNCTION_MONO_NFCHOA_25D(x0,xs,src,f,conf) returns the
%   driving signal for the given secondary source and desired source type (src).
%   The driving signal is calculated for the NFC-HOA 2.5 dimensional case in the
%   frequency domain.
%
%   References:
%
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wave_field_mono_nfchoa_25d,
%             driving_function_imp_nfchoa_25d

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargposition(xs);
xs = position_vector(xs);
isargpositivescalar(f);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% phase of omega
phase = conf.phase;
% xref
xref = position_vector(conf.xref);
% speed of sound
c = conf.c;
% positions and angles of secondary sources
x0 = x0(:,1:3);
al=atan2(x0(:,2),x0(:,1));
% number of secondary sources
M=length(x0);
% max order of spherical harmonics
if(isodd(M))
    N=(M-1)/2;
else
    N=floor((M+1)/2);
end
% source angle
al_pw=atan2(xs(2),xs(1));
% radius of array
R = norm(xref-x0(1,1:3));


%% ===== Computation ====================================================
% Calculate the driving function in frequency domain

D=zeros(1,M);
% wavenumber
k = 2*pi*f/c;


if strcmp('pw',src)

    % ===== PLANE WAVE ===============================================
    for m=0:M-1
        for n = -N : N
            D(m+1) = D(m+1) + 4.*pi .* ...
                ( 1i.^(-abs(n)) .* exp(-1i.*n.*al_pw) ) ./ ...
                ( -1i .* k .* sphbesselh(abs(n),2,k.*R) ) .*exp(1i.*n.*al(m+1));
        end
    end

    elseif strcmp('ps',src)

    % ===== POINT SOURCE =============================================
    to_be_implemented(mfilename);

    elseif strcmp('fs',src)

    % ===== FOCUSED SOURCE ===========================================
    to_be_implemented(mfilename);

    else
    % No such source type for the driving function
    error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
 end

