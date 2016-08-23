function [d,x0,xv] = driving_function_imp_localwfs(x0,xs,src,conf)
%DRIVING_FUNCTION_IMP_LOCALWFS returns the driving signal d for local WFS
%
%   Usage: [d,x0,xv] = driving_function_mono_localwfs(x0,xs,src,conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [nx6]
%       xs          - position of virtual source or direction of plane
%                     wave / m [1x3]
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'ls' - line source
%                         'fs' - focused source
%
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       d           - driving function signal [nx1]
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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargxs(xs);
isargchar(src);
isargstruct(conf);


%% ===== Configuration ==================================================
virtualconf = conf;
virtualconf.secondary_sources.size = conf.localsfs.vss.size;
virtualconf.secondary_sources.center = conf.localsfs.vss.center;
virtualconf.secondary_sources.geometry = conf.localsfs.vss.geometry;
virtualconf.secondary_sources.number = conf.localsfs.vss.number;
virtualconf.usetapwin = conf.localsfs.usetapwin;
virtualconf.tapwinlen = conf.localsfs.tapwinlen;
virtualconf.wfs = conf.localsfs.wfs;
method = conf.localsfs.method;

N = conf.N;
fs = conf.fs;


%% ===== Computation ====================================================
if strcmp('fs',src)
  error(['%s: %s is not a supported method source type! Try to use a point', ...
    ' source, if the source is inside the secondary source array but not', ...
    ' inside the virtual secondary source array'], upper(mfilename),src);
end

% Determine driving functions of virtual array with different sfs methods
switch method
  case 'wfs'
    % === Use WFS for the virtual secondary sources ===
    % Create virtual source array
    xv = virtual_secondary_source_positions(x0,xs,src,conf);
    % selection of virtual secondary source
    xv = secondary_source_selection(xv, xs, src);
    % Optional tapering
    xv = secondary_source_tapering(xv,virtualconf);
    
    % Calculate pre-equalization filter if required ( hpre(t) * hpre(-t) )
    pulse = wfs_preequalization(dirac_imp(),conf);    
    pulse = conv( wfs_preequalization(dirac_imp(),virtualconf), ...
      pulse(end:-1:1) );
    
    % Source position
    xs = repmat(xs(1:3),[size(xv,1) 1]);
    
    % Get the delay and weighting factors for the virtual secondary sources
    if strcmp('pw',src)
      % === Plane wave =====================================================
      % Direction of plane wave
      nk = bsxfun(@rdivide,xs,vector_norm(xs,2));
      % Delay and amplitude weight
      [tauv,wv] = driving_function_imp_wfs_pw(xv(:,1:3),xv(:,4:6),nk,conf);
      
    elseif strcmp('ps',src)
      % === Point source ===================================================
      % Delay and amplitude weight
      [tauv,wv] = driving_function_imp_wfs_ps(xv(:,1:3),xv(:,4:6),xs,conf);
      
    elseif strcmp('ls',src)
      % === Line source ====================================================
      % Delay and amplitude weight
      [tauv,wv] = driving_function_imp_wfs_ls(xv(:,1:3),xv(:,4:6),xs,conf);
      
    elseif strcmp('fs',src)
      % === Focused source =================================================
      % Delay and amplitude weight
      [tauv,wv] = driving_function_imp_wfs_fs(xv(:,1:3),xv(:,4:6),xs,conf);
    else
      error('%s: %s is not a known source type.',upper(mfilename),src);
    end
    
    % Neglect virtual secondary sources with zero weight
    selector = wv~=0 & xv(:,7)~=0;
    xv = xv( selector, :);
    tauv = tauv( selector );
    wv = wv( selector );
    
    % Select secondary sources based on the positions of the virtual
    % secondary sources
    x0 = secondary_source_selection(x0, xv(:,1:6), 'vss');
    
    % Initialize
    Nv = size(xv,1);
    N0 = size(x0,1);
    
    tau0 = inf(N0, Nv);  % rows for ss; columns for vss
    w0 = zeros(N0, Nv);  % rows for ss; columns for vss
    
    % interate over all virtual secondary sources
    idx = 1;
    for xvi = xv'
      % Select active source for one focused source
      [x0s, xdx] = secondary_source_selection(x0,xvi(1:6)','fs');
      if ~isempty(x0s) && xvi(7) > 0
        % Focused source position
        xs = repmat(xvi(1:3)',[size(x0s,1) 1]);
        % Delay and weights for single focused source
        [tau0(xdx,idx),w0(xdx,idx)] = driving_function_imp_wfs_fs( ...
          x0s(:,1:3),x0s(:,4:6),xs,conf);
        % Optional tapering
        x0s = secondary_source_tapering(x0s,conf);
        % Apply secondary sources' tapering and possibly virtual secondary
        % sources' tapering to weighting matrix
        w0(xdx,idx) = w0(xdx,idx).*wv(idx).*x0s(:,7).*xvi(7);
        % add up delay of secondary sources and virtual secondary sources
        tau0(xdx,idx) = tau0(xdx,idx) + tauv(idx);
      end
      idx = idx + 1;
    end    
    
    % SWITCH DIMENSIONS OF WEIGHTS AND DELAYS
    w0 = w0.';      % NOW: rows for vss; columns for ss
    tau0 = tau0.';  % NOW: rows for vss; columns for ss    
    % select non-zero weights
    selector = w0 ~= 0;
    % Remove delay offset, in order to begin always at t=0 with the first wave 
    % front at any secondary source
    tau0 = tau0 - min(tau0(selector));     
    % Shift and weight prototype driving function
    pulse = repmat( [pulse; zeros(N-length(pulse),1)] , 1, sum(selector(:)) );
    pulse = delayline(pulse, tau0(selector)*fs, w0(selector), conf);
    % Compose impulse responses
    d = zeros(N, N0);
    kdx = 1;
    for idx=1:N0
      l = sum( w0(:, idx) ~= 0 );
      if l > 0
        % sum up prototypes which belong to the idx'th secondary source
        d(:, idx) = sum( pulse(:, kdx:kdx+l-1), 2 );
        kdx = kdx + l;
      end
    end
  otherwise
    error('%s: %s is not a supported method for time domain localsfs!', ...
      upper(mfilename),method);
end

