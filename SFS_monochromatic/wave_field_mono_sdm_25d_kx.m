function [x,y,P] = wave_field_mono_sdm_25d_kx(X,Y,xs,src,f,L,conf)
%WAVE_FIELD_SDM_WFS_25D_KX simulates the wave field of a given source for 25D SDM
%IN THE SPATIAL FREQUENCY DOMAIN
%   Usage: [x,y,P] = wave_field_mono_sdm_25d_kx(X,Y,xs,src,f,L,[conf])
%
%   Input parameters:
%       X           - [xmin,xmax]
%       Y           - [ymin,ymax]
%       xs          - position of point source (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       f           - monochromatic frequency (Hz)
%       L           - array length (m)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x           - corresponding x axis
%       y           - corresponding y axis
%       P           - Simulated wave field
%
%   WAVE_FIELD_MONO_SDM_25D_KX(X,Y,xs,src,f,L,conf) simulates a wave field of
%   the given source type (src) using a SDM 2.5 dimensional driving function
%   in the spectro-temporal freqeuncy domain. 
%   To plot the result use plot_wavefield(x,y,P).
%
%   NOTE: due to numerical problems with the fft and the bessel functions needed
%   in SDM (which resulted in an imaginary part which is hundreds of orders
%   greater/smaller than the real part) the FFT is done by hand in this
%   function. This results in a longer time to run this function. If you haven't
%   that time and you can try the large argument approximation of the
%   bessel functions, which will result in a wrong evanescent part of the wave
%   field.
%
%   References:
%       Spors2010 - Reproduction of Focused Sources by the Spectral Division
%           Method
%       Spors2010 - Analysis and Improvement of Pre-equalization in
%       2.5-Dimensional Wave Field Synthesis
%
%   see also: plot_wavefield, wave_field_mono_wfs_25d

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
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(X,Y);
xs = position_vector(xs);
isargpositivescalar(L,f);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Array position (m)
X0 = conf.X0;
% Loudspeaker distcane
% NOTE: if dLS <= dx, than we have a continous loudspeaker
dx = conf.dx0;
% Method to calculate driving function (only for non-aliased part)
withev = conf.withev;  % with evanescent waves
% Reference position for the amplitude (correct reproduction of amplitude
% at y = yref).
xref = conf.xref;
c = conf.c;
% Plotting
useplot = conf.useplot;
% xy resolution
xysamples = conf.xysamples;


%% ===== Variables ======================================================
% General
omega = 2*pi*f;
% init variables
kxrep=2*pi/dx;
Nrep=6;             % number of spectral repetitions
% Aliasing condition
kxal = omega/c;
% Factor by which kx is extended of kx = omega/c criteria
Nkx=1.5;
%kx = linspace(-Nkx*kxal,Nkx*kxal,Nkx*2000);
kx = linspace(-Nkx*kxal,Nkx*kxal,Nkx*xysamples*10);
x  = linspace(X(1),X(2),xysamples);
y  = linspace(Y(1),Y(2),xysamples);
% Indexes for evanescent contributions and propagating part of the wave field
idxpr = (( abs(kx) <= (omega/c) ));
idxev = (( abs(kx) > (omega/c) ));


%% ===== Wave field in the spectro-temporal domain =======================
%
% ========================================================================
% Secondary source model
Gkx = zeros(length(kx),length(y));
% Green's function for a point source in the spectro-temporal domain (see
% Spors2010)
%                                  ____________
%                 / -i/4 H0^(2)( \|(w/c)^2-kx^2 y )
% G_3D(kx,y,w) = <                ____________
%                 \ 1/(2pi) K0( \|kx^2-(w/c)^2 y )
%
[K,Y] = meshgrid(kx(idxpr),abs(y-X0(2)));
Gkx(idxpr,:) = -1j/4 .* besselh(0,2,sqrt( (omega/c)^2 - K.^2 ).* Y)';
if(withev)
    [K,Y] = meshgrid(kx(idxev),abs(y-X0(2)));
    Gkx(idxev,:) = 1/(2*pi) .* besselk(0,sqrt( K.^2 - (omega/c)^2).* Y)';
end

% ========================================================================
% Driving function
Dkx = driving_function_mono_sdm_25d_kx(kx,xs,src,f,conf);
% Convolution with a window representing the length L of the loudspeaker array
% FIXME: this doesn't work with evanescent waves at the moment
%w = L * sin(kx*L/2)./(kx*L/2);
%Dkx = conv2(Dkx,w,'same');


%% =======================================================================
% Reproduced field
% Pkx = Dkx * Gkx
Pkx = repmat(Dkx',1,length(y)) .* Gkx;


%% ===== Inverse spatial Fourier transformation =========================
% 
%            /
% P(x,y,w) = | Pkx(kx,y,w) * e^(-i kx x) dkx
%            /
%
P = zeros(length(y),length(x));
for n=1:length(x)
    % The following loop can be done faster by using the line below with repmat
    %for m=1:length(y)
    %    P(m,n) = sum ( Pkx(:,m) .* exp(-1j*kx*x(n))' )';
    %end
    P(:,n) = sum ( Pkx .* repmat(exp(-1j*kx*x(n))',1,xysamples),1 )';
end

% === Scale signal (at [xref yref]) ===
P = norm_wave_field(P,x,y,conf);

%% ===== Plotting ========================================================
if(useplot)
    plot_wavefield(x,y,P,conf);
end
