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

%% ===== Computation =====================================================
% Generate stimulus pusle
pulse = dirac_imp();
% Radius of array
R = norm(x0(1,1:3)-X0);
% Ambisonics order
if isempty(conf.nfchoa.order)
    % Get maximum order of spherical harmonics
    order = nfchoa_order(nls,conf);
else
    order = conf.nfchoa.order;
end

% Correct position of source for off-center arrays
xs(1:3) = xs(1:3)-X0;
[theta_src, r_src] = cart2pol(xs(1),xs(2));

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

% -------------------------------------------------------------------------
%                      ___
%                1     \
% D(phi0,w) = -------  /__     Hm(w) e^(im(phi0-phi_src))
%             2 pi r0  m=-M..M
%
% see Spors et al. (2011), eq.(4) and (5)
%--------------------------------------------------------------------------

% apply modal weighting
wm = modal_weighting(order,conf);  % modal weighting function
Hm = bsxfun(@times,Hm,wm);

% apply phase shift due to source angle
m = 0:order;
dm = bsxfun(@times,Hm,exp(-1i*m*theta_src));

% Spatial IFFT
dm = [conj(dm(:,end:-1:2)), dm];  % append coefficients for negative m
d = inverse_cht(dm,nls);
d = 1/(2*pi*R)*real(d);
