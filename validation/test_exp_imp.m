function status = test_exp_imp(modus)
%TEST_EXP_IMP tests the different expansions of sound fields in time-domain
%
%   Usage: status = test_exp_imp(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false
%
%   TEST_EXP_IMP(modus) checks, if the circular basis expansions for plane waves
%   and point sources are working. The circular basis expansions are converted
%   to plane wave decompositions. Additionally, modal weighting functions are
%   tested. Optionally, sound field plots of the plane wave decompositions are
%   used for verification.

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


status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
% Parameters
conf = SFS_config;
conf.plot.loudspeakers = 0;  % do not plot loudspeakers
conf.plot.usedb = 1;  % use decibel for sound field plot
conf.modal_window_parameter = 2.0;  % parameter for kaiser window

c = conf.c;  % speed of sound

% for sound field plots
X = [-2, 2];
Y = [-2, 2];
Z = 0;

fhp = 1000;  % high-pass filter corner-frequency for pwd of point source

% plane waves with equi-angular distribution
Npw = 1024;
phi0 = (0:Npw-1).'*2*pi/Npw;
x0 = [cos(phi0) sin(phi0)];
x0(:,3) = 0;
x0(:,4:6) = x0(:,1:3);
x0(:,7) = 1;

% test scenarios
scenarios = { ...
  'pw', [ 0.0 -1.0 0.0],  5, 'rect'   , [0.0 0.0 0.0]
  'pw', [ 0.0 -1.0 0.0], 15, 'rect'   , [0.0 0.0 0.0]
  'pw', [ 0.0 -1.0 0.0], 15, 'kaiser' , [0.0 0.0 0.0]
  'pw', [ 0.0 -1.0 0.0],  5, 'rect'   , [0.5 1.0 0.0]
  'pw', [ 0.0 -1.0 0.0], 15, 'rect'   , [0.5 1.0 0.0]
  'pw', [ 0.0 -1.0 0.0], 15, 'kaiser' , [0.5 1.0 0.0]
  'ps', [ 0.0  2.5 0.0],  5, 'rect'   , [0.0 0.0 0.0]
  'ps', [ 0.0  2.5 0.0], 15, 'rect'   , [0.0 0.0 0.0]
  'ps', [ 0.0  2.5 0.0], 15, 'kaiser' , [0.0 0.0 0.0]
  'ps', [ 0.0  2.5 0.0],  5, 'rect'   , [0.5 1.0 0.0]
  'ps', [ 0.0  2.5 0.0], 15, 'rect'   , [0.5 1.0 0.0]
  'ps', [ 0.0  2.5 0.0], 15, 'kaiser' , [0.5 1.0 0.0]
  };

%% ===== Main ============================================================

for ii=1:size(scenarios)

    src = scenarios{ii,1};  % source type
    xs = scenarios{ii,2};  % source position / direction of plane wave
    Nce = scenarios{ii,3};  % modal order
    conf.modal_window = scenarios{ii,4};  % type of modal weighting function
    xq = scenarios{ii,5};  % expansion centre

    % circular expansion coefficients
    switch src
    case 'pw'
        [pm, delay_offset] = circexp_imp_pw(xs,Nce,xq,conf);
        t = xq*xs.'./c;
    case 'ps'
        [pm, delay_offset] = circexp_imp_ps(xs,Nce,xq,fhp,conf);
        t = norm(xq-xs)./c;
    end

    % modal weighting of coefficients
    wm = modal_weighting(Nce, conf);
    pm = bsxfun(@times, pm, wm);

    % conversion to plane wave decomposition
    ppwd = pwd_imp_circexp(pm, Npw);

    if modus
      % delay plane waves according to expansion centre
      delays = -xq*x0(:,1:3).'./c;
      delay_offset = delay_offset - min(delays);
      delays = delays - min(delays);
      ppwd = delayline(ppwd,delays,1,conf);

      % sound field plot
      sound_field_imp(X,Y,Z,x0,'pw',ppwd,t+delay_offset,conf);

      % title string
      str = 'Plane wave decompostion of modally bandlimited';
      switch src
      case 'pw'
          str = sprintf('%s plane wave', str);
      case 'ps'
          str = sprintf('%s point source', str);
      end
      str = sprintf(['%s ([%1.1f %1.1f %1.1f]):\n%s-window (M=%d), ' ...
        'center of modal expansion at [%1.1f %1.1f %1.1f]'], ...
        str, xs, conf.modal_window, Nce, xq);
      title(str);

    end
end

status = true;
