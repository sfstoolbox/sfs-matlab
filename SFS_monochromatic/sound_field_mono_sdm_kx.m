function varargout = sound_field_mono_sdm_kx(X,Y,Z,xs,src,f,conf)
%SOUND_FIELD_MONO_SDM_KX simulates the sound field of a given source for SDM
%in the kx domain
%
%   Usage: [P,x,y,z] = sound_field_mono_sdm_kx(X,Y,Z,xs,src,f,conf)
%
%   Input parameters:
%       X           - x-axis / m; [xmin,xmax]
%       Y           - y-axis / m; [ymin,ymax]
%       Z           - z-axis / m; single value
%       xs          - position of point source / m
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       f           - monochromatic frequency / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - Simulated sound field
%       x           - corresponding x axis / m
%       y           - corresponding y axis / m
%       z           - corresponding z axis / m
%
%   SOUND_FIELD_MONO_SDM_KX(X,Y,Z,xs,src,f,conf) simulates a monochromatic sound
%   field of the given source type (src) synthesized with the spectral devision
%   method (SDM). Note, that the linaer secondary sources are placed automatically
%   on a line parrallel to the x-axis accordingly to conf.secondary_sources.center.
%   The field can only be calculated in the xy-plane, meaning only Z=0 is allowed.
%
%   To plot the result use:
%   plot_sound_field(P,X,Y,Z,conf);
%   or simple call the function without output argument:
%   sound_field_mono_sdm_kx(X,Y,Z,xs,src,f,conf)
%%
%   NOTE: due to numerical problems with the fft and the bessel functions needed
%   in SDM (which resulted in an imaginary part which is hundreds of orders
%   greater/smaller than the real part) the FFT is done by hand in this
%   function. This results in a longer time to run this function. If you haven't
%   that time and you can try the large argument approximation of the
%   bessel functions, which will result in a wrong evanescent part of the sound
%   field.
%
%   References:
%       S. Spors, J. Ahrens (2010) - "Analysis and Improvement of
%       Pre-equalization in 2.5-Dimensional Wave Field Synthesis", 128th AES
%       Conv.
%
%   See also: plot_sound_field, sound_field_mono_sdm

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
nargmin = 7;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z);
if length(Z)>1 || Z~=0
    error('%s: SDM in the kx domain works only for Z=0 at the moment.', ...
        upper(mfilename));
end
isargxs(xs);
isargpositivescalar(f);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ==================================================
% Array position (m)
X0 = conf.secondary_sources.center;
withev = conf.sdm.withev;  % with evanescent waves
c = conf.c;
% xy resolution
resolution = conf.resolution;
useplot = conf.plot.useplot;


%% ===== Variables ======================================================
% General
omega = 2*pi*f;
% Aliasing condition
kxal = omega/c;
% Factor by which kx is extended of kx = omega/c criteria
Nkx=1.5;
kx = linspace(-Nkx*kxal,Nkx*kxal,Nkx*resolution*10);
% Create axes
x = linspace(X(1),X(2),resolution);
y = linspace(Y(1),Y(2),resolution);
z = Z;
% Indexes for evanescent contributions and propagating part of the sound field
idxpr = (( abs(kx) <= (omega/c) ));
idxev = (( abs(kx) > (omega/c) ));


%% ===== Wave field in the spectro-temporal domain =======================
%
% === Secondary source model ===
Gkx = zeros(length(kx),length(y));
% Green's function for a point source in the spectro-temporal domain, see
% Spors (2010)
%                                  ____________
%                 / -i/4 H0^(2)( \|(w/c)^2-kx^2 y )
% G_3D(kx,y,w) = <                ____________
%                 \ 1/(2pi) K0( \|kx^2-(w/c)^2 y )
%
[kk,yy] = meshgrid(kx(idxpr),abs(y-X0(2)));
Gkx(idxpr,:) = -1j/4 .* besselh(0,2,sqrt( (omega/c)^2 - kk.^2 ).* yy).';
if(withev)
    [kk,yy] = meshgrid(kx(idxev),abs(y-X0(2)));
    Gkx(idxev,:) = 1/(2*pi) .* besselk(0,sqrt( kk.^2 - (omega/c)^2).* yy).';
end

% ========================================================================
% Driving function
Dkx = driving_function_mono_sdm_kx(kx,xs,src,f,conf);
% Convolution with a window representing the length L of the loudspeaker array
% FIXME: this doesn't work with evanescent waves at the moment
%w = L * sin(kx*L/2)./(kx*L/2);
%Dkx = conv2(Dkx,w,'same');


%% =======================================================================
% Reproduced field
% Pkx = Dkx * Gkx
Pkx = repmat(Dkx.',1,length(y)) .* Gkx;


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
    P(:,n) = sum ( Pkx .* repmat(exp(-1j*kx*x(n))',1,resolution),1 )';
end

% return parameter
if nargout>0, varargout{1}=P; end
if nargout>1, varargout{2}=x; end
if nargout>2, varargout{3}=y; end
if nargout>3, varargout{4}=z; end


%% ===== Plotting ========================================================
if nargout==0 || useplot
    plot_sound_field(P,X,Y,Z,conf);
end
