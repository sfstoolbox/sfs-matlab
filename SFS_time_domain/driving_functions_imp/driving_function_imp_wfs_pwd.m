function [d,delay_offset] = driving_function_imp_wfs_pwd(x0,ppwd,xq,conf)
%DRIVING_FUNCTION_IMP_WFS_PWD returns the WFS driving signal for a plane wave 
%expansion
%
%   Usage: [d,delay_offset] = driving_function_imp_wfs_pwd(x0,ppwd,[xq],conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [N0x6]
%       ppwd        - plane wave coefficients [N x Npw]
%       xq          - centre of plane wave expansionl, default = [0,0,0];
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       d           - driving function signals [nxN0]
%       
%       
%       x0          - position, direction, and weights of the real secondary
%                     sources / m [nx7]
%       xv          - position, direction, and weights of the virtual secondary
%                     sources / m [mx7]
%
%   References:
%       S. Spors (2010) - "Local Sound Field Synthesis by Virtual Secondary
%                          Sources", 40th AES
%
%   See also: plot_sound_field, sound_field_mono_wfs

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargmatrix(ppwd);
if nargin == nargmin
    conf = xq;
    xq = [0, 0, 0];
else
    isargxs(xq);
end
isargstruct(conf);


%% ===== Computation ====================================================
% some sizes
Npw = size(ppwd, 2);
N0 = size(x0,1);

% apply integrations weights to pwd
ppwd = ppwd./Npw;

% Calculate pre-equalization filter if required hpre(t)
[hpre,delay_offset] = wfs_preequalization(dirac_imp(),conf);

% apply hpre to pwd if its more efficient (less plane wave than loudspeakers)
if Npw <= N0
    ppwd = convolution(ppwd,hpre);
end
N = size(ppwd,1);

% shift coordinates to expansion center
x0(:,1:3) = bsxfun(@minus,x0(:,1:3),xq);
conf.xref = [0 0 0];

% initialise delay and weights
tau0 = inf(Npw,N0);  % rows for plane waves; columns for ss
w0 = zeros(Npw,N0);  % rows for plane waves; columns for ss

% interate over all plane wave directions
phipw = (0:Npw-1)*2*pi/Npw;
npw = [cos(phipw); sin(phipw)];  % [2 x Npw]
npw(3,:) = 0;
ndx = 0;
for nk = npw
    ndx = ndx + 1;
    % Select active source for one focused source
    [x0s,xdx] = secondary_source_selection(x0,nk.','pw');

    if ~isempty(x0s)
        % Focused source position
        xs = repmat(nk.',[size(x0s,1) 1]);
        % Delay and weights for single plane wave
        [tau0(ndx,xdx),w0(ndx,xdx)] = driving_function_imp_wfs_pw( ...
            x0s(:,1:3),x0s(:,4:6),xs,conf);
        % Optional tapering
        x0tmp = secondary_source_tapering(x0s,conf);
        wtap = x0tmp(:,7)./x0s(:,7);  % x0
        % Apply secondary sources' tapering to weighting matrix
        w0(ndx,xdx) = w0(ndx,xdx).*wtap.';
    end
end
% Remove delay offset, in order to begin always at t=0 with the first wave 
% front at any secondary source
delay_offset = delay_offset - min(tau0(w0~=0));
tau0 = tau0 - min(tau0(w0~=0)); 

% Compose impulse responses
d = zeros(N,N0);
for ndx = 1:Npw  
    % select only secondary source with non-zero weights
    xsel = w0(ndx,:) ~= 0;
    
    [tmp,delayline_delay] = ...
        delayline(ppwd(:,ndx),tau0(ndx,xsel),w0(ndx,xsel),conf);
    d(:,xsel) = d(:,xsel) + tmp;
end
% add delay of delayline
delay_offset = delay_offset + delayline_delay;

% apply hpre to driving signals if its more efficient
if Npw > N0
    d = convolution(d,hpre);
end
d = d(1:conf.N, :);
