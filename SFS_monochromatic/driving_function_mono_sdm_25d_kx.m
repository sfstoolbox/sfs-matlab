function [Dkx] = driving_function_mono_sdm_25d_kx(kx,xs,src,f,conf)
%DRIVING_FUNCTION_MONO_SDM_25D_KX returns the driving signal D for 2.5D SDM
%IN THE SPATIAL FREQUENCY (kx-)DOMAIN
%
%   Usage: Dkx = driving_function_mono_sdm_25d_kx(kx,xs,src,f,[conf])
%
%   Input parameters:
%       kx          - position of the secondary source in the spectro-temporal
%                     domain (m)
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
%       Dkx         - driving function signal (n x n)
%
%   DRIVING_FUNCTION_MONO_SDM_25D_KX(x0,xs,f,src,conf) returns the
%   driving signal for the given secondary source and desired source type (src).
%   The driving signal is calculated for the SDM 2.5 dimensional case in the
%   spatial frequency domain.

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
isargmatrix(kx);
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
phase = conf.phase;
xref = position_vector(conf.xref);
X0 = conf.X0;
c = conf.c;
withev = conf.withev;


%% ===== Variables =======================================================
% Omega
omega = 2*pi*f;
% Indexes for evanescent contributions and propagating part of the wave field
idxpr = (( abs(kx) <= (omega/c) ));
idxev = (( abs(kx) > (omega/c) ));


%% ===== Spectrum of driving function ====================================
Dkx = zeros(1,length(kx));

if strcmp('pw',src)

    % ===== PLANE WAVE ===================================================
    al_pw=atan2(xs(2),xs(1));
    kpwx=omega/c*cos(al_pw);
    kpwy=omega/c*sin(al_pw);
    idx = find(kx>=kpwx,1,'first');
    
    Dkx(idx) = 4i*exp(-1i*kpwy*xref(2))./besselh(0,2,kpwy*xref(2));
    
    % spectral repetitions
    %     Nrep=6;
    %     kxrep=2*pi/0.2;
    %     if(Nrep>0)
    %         for n=-Nrep:Nrep
    %             if(n~=0)
    %                 kxp=kpwx-n*kxrep;
    %     
    %                 idx=find(kx>=kxp,1,'first');
    %                     Dkx(idx) = 4i*exp(-1i*kpwy*xref(2))./besselh(0,2,kpwy*xref(2));
    %             end
    %         end
    %     end
    

elseif strcmp('ps',src)

    % ===== POINT SOURCE =================================================
    to_be_implemented(mfilename);

elseif strcmp('fs',src)

    % ===== FOCUSED SOURCE ===============================================
    %
    % Exact driving function for 25D synthesis
    %
    % D_25D(kx,w) = e^(i kx xs) ...
    %                                   ____________
    %                         H0^(2)( \|(w/c)^2-kx^2 |yref-ys| )
    %                     / - --------------_-_-_-_-_-_---------, |kx|<|w/c|
    %                     |      H0^(2)( \|(w/c)^2-kx^2 yref ) 
    %                    <        ____________
    %                     | K0( \|kx^2-(w/c)^2 |yref-ys| )
    %                     \ ----------_-_-_-_-_-_---------,       |kx|>|w/c|
    %                          K0( \|kx^2-(w/c)^2 yref )
    %
    % NOTE: the phase term e^(i phase) is only there in order to be able
    %       to simulate different time steps
    %
    Dkx(idxpr) =  exp(1i*kx(idxpr)*xs(1)) .* ...
        besselh(0,2,sqrt( (omega/c)^2 - kx(idxpr).^2 )*abs(xref(2)-xs(2))) ./ ...
        besselh(0,2,sqrt( (omega/c)^2 - kx(idxpr).^2 )*abs(xref(2)-X0(2))) .* ...
        exp(1i*phase);
    if(withev)
        Dkx(idxev) =  exp(1i*kx(idxev)*xs(1)) .* ...
            besselk(0,sqrt(kx(idxev).^2 - (omega/c).^2)*abs(xref(2)-xs(2))) ./ ...
            besselk(0,sqrt(kx(idxev).^2 - (omega/c).^2)*abs(xref(2)-X0(2))) .* ...
            exp(1i*phase);
    end

else
    % No such source type for the driving function
    error('%s: src has to be one of "pw", "ps", "fs"!',upper(mfilename));
end
