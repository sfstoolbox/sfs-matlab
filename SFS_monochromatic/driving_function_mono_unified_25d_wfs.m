function [D, xPCS] = driving_function_mono_unified_25d_wfs(x0,xs,dx0,src,f,conf)
%DRIVING_FUNCTION_MONO_UNIFIED_25D_WFS returns the driving signal D for WFS
%
%   derived driving_function_mono_wfs(x0,xs,src,f,conf)
%
%   Usage: D = driving_function_mono_unified_25d_wfs(x0,xs,dx0,src,f,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx6]
%       xs          - position of virtual source or direction of plane
%                     wave / m [1x3] or [1x6]
%       dx0         - amplitude factor of unified 2.5D WFS framework
%                     to obtain amplitude correct synthesis at desired
%                     locations, i.e. along a definable reference curve,
%                     this is primary source AND x0 dependent,
%                     one dx0 per one x0, thus [nx1]
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'ls' - line source
%                         'fs' - focused source
%       f           - temporal frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function temporal spectrum [nx1]
%       xPCS        - locations of amplitude correct synthesis [nx3]
%
%   References:
%   [Sch17] Frank Schultz, Gergely Firtha, Peter Fiala, Sascha Spors (2017):
%   "Wave Field Synthesis Driving Functions for Large-Scale Sound
%   Reinforcement Using Line Source Arrays." In: Proc. of 142nd Audio Eng.
%   Soc. Conv., Berlin, #9722.
%
%   [Fir17] Gergely Firtha, Peter Fiala, Frank Schultz, Sascha Spors (2017):
%   "Improved Referencing Schemes for 2.5D Wave Field Synthesis Driving
%   Functions." In: IEEE/ACM Trans. Audio Speech Language Process.,
%   DOI 10.1109/TASLP.2017.2689245
%
%   DRIVING_FUNCTION_MONO_UNIFIED_25D_WFS(x0,xs,dx0,src,f,conf)
%   returns the driving signal for
%   the given secondary sources and desired source type (src) for the
%   given temporal frequency using the unified 2.5D WFS framework
%
%   See also: plot_sound_field, sound_field_mono_wfs_25d,
%             driving_function_imp_wfs_25d

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
nargmin = 6;
nargmax = 6;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargxs(xs);

%isargdx0(dx0); TBD!!!

isargpositivescalar(f);
isargchar(src);
isargstruct(conf);


%% ===== Computation ====================================================
% Calculate the driving function in temporal-frequency domain

% Secondary source positions and directions
nx0 = x0(:,4:6);
x0 = x0(:,1:3);

% Source position/direction/orientation
xs = repmat(xs,[size(x0,1) 1]);

%%
% unified 2.5D WFS framework driving function, cf. [(1), Sch17], [(47), Fir17]
% get driving signals
omega = 2*pi*f;
c = conf.c;
w_c = omega/c;
pre = -sqrt(8*pi/(1i*w_c));

if strcmp('pw',src)
    % === Plane wave ===
    % directional derivative [(3)&(4), Sch17]:
    nk = bsxfun(@rdivide,xs,vector_norm(xs(:,1:3),2)); % propagating direction of plane wave
    dP_dn = (-1i*w_c) .*...
        vector_product(nk,nx0,2) .*...
        exp(-1i.*w_c.*vector_product(nk,x0,2));
    
    xPCS = x0 + nk.*repmat(dx0,[1,3]); %get the locations/positions of
    %amplitude correct synthesis, [(8), Sch17]
    
elseif strcmp('ps',src)
    % === Point source ===
    % directional derivative [(15)&(16), Sch17]
    r = vector_norm(x0-xs,2); % r = |x0-xs|
    dP_dn = (-1i*w_c) .*...
        (vector_product(x0-xs,nx0,2)./r) .*...
        (exp(-1i.*w_c.*r)./(4*pi*r));

    xPCS = x0 + (x0-xs) .* repmat(dx0./(r-dx0),[1,3]); %get the locations
    %/positions of amplitude correct synthesis, [(20), Sch17],
    %[(33)&(34), Fir17]
    
elseif strcmp('ls',src)
    % === Line source ===
    % directional derivative [(10), Sch17]
    r = vector_norm(x0-xs,2); % r = |x0-xs|    
    dP_dn = (1i*w_c)/4 .*...
        (vector_product(x0-xs,nx0,2)./r) .*...
        besselh(1,2,w_c.*r);

    xPCS = x0 + ((x0-xs)./r).*repmat(dx0,[1,3]); %get the locations/
    %positions of amplitude correct synthesis, [(14), Sch17] 

elseif strcmp('fs',src)
    % === Focused source ===
    error('%s: %s is not implemented yet.',upper(mfilename),src);
else
    error('%s: %s is not a known source type.',upper(mfilename),src);
end
% put all together for driving function
D = pre .*  sqrt(dx0) .* dP_dn;

end
