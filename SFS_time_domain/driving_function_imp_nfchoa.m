function [d,dm,delay_offset] = driving_function_imp_nfchoa(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_NFCHOA calculates the NFC-HOA driving function
%
%   Usage: [d,dm,delay_offset] = driving_function_imp_nfchoa(x0,xs,src,conf)
%
%   Input parameters:
%       x0      - position  and direction of secondary sources [N0x7] / m
%       xs      - position of virtual source or direction of plane wave / m
%       src     - source type of the virtual source
%                     'pw' - plane wave (xs, ys are the direction of the
%                            plane wave in this case)
%                     'ps' - point source
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       d            - matrix of driving signals [NxN0]
%       dm           - matrix of driving funtion in spherical/circular domain
%                      [NxM]
%       delay_offset - delay add by driving function / s
%
%   DRIVING_FUNCTION_IMP_NFCHOA(x0,xs,src,conf) returns the
%   driving function of NFC-HOA for the given source type and position,
%   and loudspeaker positions.
%
%   References:
%     Spors, S., Kuscher, V., Ahrens, J. (2011) - "Efficient realization of
%         model-based rendering for 2.5-dimensional near-field compensated
%         higher order Ambisonics", IEEE Workshop on Applications of Signal
%         Processing to Audio and Acoustics (WASPAA), pp. 61-64,
%         https://doi.org/10.1109/ASPAA.2011.6082325
%     Schultz, F. and Spors, S. (2014) - "Comparing Approaches to the Spherical
%         and Planar Single Layer Potentials for Interior Sound Field
%         Synthesis," Acta Acustica united with Acustica, pp. 900-911,
%         https://doi.org/10.3813/AAA.918769
%     Gumerov, N. and Duraiswami, R. (2004) - "Fast multipole methods for the 
%         Helmholtz equation in three dimensions," Elsevier, Oxford
%
%   See also: driving_function_imp_nfchoa_ps, sound_field_imp_nfchoa

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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0)
isargxs(xs);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ==================================================
nls = size(x0,1);
N = conf.N;
X0 = conf.secondary_sources.center;
t0 = conf.t0;
c = conf.c;
dimension = conf.dimension;

%% ===== Computation =====================================================
% Generate stimulus pusle
pulse = dirac_imp();
% Correct position of loudspeakers for off-center arrays
x0(:,1:3) = bsxfun(@minus, x0(:,1:3), X0);
% Radius of array
R = norm(x0(1,1:3));
% Ambisonics order
if isempty(conf.nfchoa.order)
    % Get maximum order of spherical harmonics
    order = nfchoa_order(nls,conf);
else
    order = conf.nfchoa.order;
end

% Correct position of source for off-center arrays
xs(1:3) = xs(1:3)-X0;
[phi_src, ~, r_src] = cart2sph(xs(1),xs(2),xs(3));

% Compute impulse responses of modal filters
Hm = [repmat(pulse,[1 order+1]); zeros(N-length(pulse),order+1)];
for n=1:order+1

    % Get the second-order sections for the different virtual sources
    if strcmp('pw',src)
        % === Plane wave =================================================
        [sos,g] = driving_function_imp_nfchoa_pw(n-1,R,conf);
    elseif strcmp('ps',src)
        % === Point source ===============================================
        [sos,g] = driving_function_imp_nfchoa_ps(n-1,R,r_src,conf);
    elseif strcmp('ls',src)
        % === Line source ================================================
        [sos,g] = driving_function_imp_nfchoa_ls(n-1,R,r_src,conf);
    elseif strcmp('fs',src)
        % === Focussed source ============================================
        [sos,g] = driving_function_imp_nfchoa_fs(n-1,R,r_src,conf);
    else
        error('%s: %s is not a known source type.',upper(mfilename),src);
    end

    % Apply them by a bilinear transform and filtering
    [b,a] = bilinear_transform(sos,conf);
    for ii=1:length(b)
        Hm(:,n) = filter(b{ii},a{ii},Hm(:,n));
    end
    Hm(:,n) = Hm(:,n)*g;  % apply gain factor
end

% Delay_offset
if strcmp('system',t0)
    delay_offset = 0;
elseif strcmp('source',t0)
    switch src
    case 'pw'
        delay_offset = R/c;
    case 'ps'
        delay_offset = (R-r_src)/c;
    end
end

% Apply modal weighting
wm = modal_weighting(order,conf);
Hm = bsxfun(@times,Hm,wm);

if strcmp('2.5D',dimension)
    % -------------------------------------------------------------------------
    %                      ___
    %                1     \
    % D(phi0,w) = -------  /__     Hm(w) e^(im(phi0-phi_src))
    %             2 pi r0  m=-M..M
    %
    % See Spors et al. (2011), eq. (4)
    %--------------------------------------------------------------------------
    
    % Apply phase shift due to source angle
    m = 0:order;
    dm = bsxfun(@times,Hm,exp(-1i*m*phi_src));
    
    % Inverse circular harmonics transform
    dm = [conj(dm(:,end:-1:2)) dm];  % append coefficients for negative m
    d = inverse_cht(dm,nls);
    d = 1/(2*pi*R)*real(d);

elseif strcmp('3D',dimension)
    % -------------------------------------------------------------------------
    %                    ___          ___
    %               1    \            \       -m               m
    % D(x0,w) =  ------  /__   Hn(w)  /__    Y  (thetaS,phiS) Y (theta0,phi0)
    %             R.^2  n=0..M      m=-n..n   n                n
    %
    % See Schultz and Spors (2014), eq. (A3) and (A6).
    % Equivalent expression, see Gumerov and Duraiswami (2004), eq. (2.1.70):
    %                    ___
    %               1    \            2n+1
    % D(x0,w) =  ------  /__   Hn(w) ------ P (cos(THETA))
    %             R.^2  n=0..M        4*pi   n
    % 
    % with P_n (cos(THETA)) being the nth-order Legendre polynomial depending
    % on the angle THETA between xs and x0. The sum corresponds to an inverse
    % Legendre transform (ILT) of Hn(w) weighted by 1/2pi.
    %--------------------------------------------------------------------------

    cosTHETA = xs*x0(:,1:3).'./r_src./R;
    dm = Hm;
    d = inverse_lt(dm,cosTHETA)./(2*pi*R.^2);  % ILT
end
