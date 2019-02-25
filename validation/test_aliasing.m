function status = test_aliasing(modus)
%TEST_ALIASING tests the prediction of aliasing frequency
%
%   Usage: status = test_driving_functions(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status - true or false

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2018 SFS Toolbox Developers                             *
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


% TODO: add mode to save data as reference data
status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
conf = SFS_config;
conf.plot.usenormalisation = false;
conf.usetapwin = false;

conf.nfchoa.order = 25;

conf.localwfs_sbl.order = 25;
conf.localwfs_sbl.Npw = 1024;
conf.localwfs_sbl.fc = 500;

conf.localwfs_vss.geometry = 'circular';
conf.localwfs_vss.size = 1.0;
conf.localwfs_vss.number = 16;
conf.localwfs_vss.usetapwin = false;
conf.localwfs_vss.consider_secondary_sources = false;
conf.localwfs_vss.consider_target_field = false;

N0gt = 256;
Nvgt = 256;
f = 1500;

% test scenarios
scenarios = { ...
    'WFS'       , 'linear',   'pw', [ 0.5 -1.0  0.0]; ...
    'WFS'       , 'linear',   'ps', [ 0.0  3.0  0.0]; ...
    'WFS'       , 'linear',   'fs', [ 0.0 -1.0  0.0  0.0 -1.0  0.0]; ...
    'LWFS-SBL'  , 'linear',   'pw', [ 0.5 -1.0  0.0]; ...
    'LWFS-SBL'  , 'linear',   'ps', [ 0.0  1.0  0.0]; ...
    'LWFS-VSS'  , 'linear',   'pw', [ 0.5 -1.0  0.0]; ...
    'LWFS-VSS'  , 'linear',   'ps', [ 0.0  1.0  0.0]; ...
    'WFS'       , 'circular', 'pw', [ 0.5  0.5  0.0]; ...
    'WFS'       , 'circular', 'ps', [ 0.0  2.5  0.0]; ...
    'WFS'       , 'circular', 'fs', [ 0.0  0.5  0.0  0.0 -1.0  0.0]; ...
    'LWFS-SBL'  , 'circular', 'pw', [ 0.5  0.5  0.0]; ...
    'LWFS-SBL'  , 'circular', 'ps', [ 0.0  2.5  0.0]; ...
    'LWFS-VSS'  , 'circular', 'pw', [ 0.5  0.5  0.0]; ...
    'LWFS-VSS'  , 'circular', 'ps', [ 0.0  2.5  0.0]; ...
    'HOA'       , 'circular', 'pw', [ 0.5  0.5  0.0]; ...
    'HOA'       , 'circular', 'ps', [ 0.0  2.5  0.0]; ...
};

% Start testing
for ii=1:size(scenarios)

    method = scenarios{ii,1};
    conf.secondary_sources.geometry = scenarios{ii,2};
    
    % set scenario
    switch scenarios{ii,2}
        case 'linear'
            X = [-1.55 1.55];
            Y = [-2.85 0.15];
            Z = 0;
            conf.xref = [0.5 -1.5 0];
            conf.secondary_sources.size = 3;
            conf.secondary_sources.number = 16;
        case {'circular', 'circle'}
            X = [-1.55 1.55];
            Y = [-1.55 1.55];
            Z = 0;
            conf.xref = [0.5, -0.75, 0];
            conf.secondary_sources.size = 3;
            conf.secondary_sources.number = 36;
    end
    xref = conf.xref;
    conf.localwfs_vss.center = conf.xref;
    
    
    % evaluation coordinates for aliasing frequency
    conf.resolution = 51;
    [xS,yS,~] = xyz_grid(X,Y,Z,conf);
    xvec = [xS(:), yS(:)];
    xvec(:,3) = 0;

    src = scenarios{ii,3};
    if strcmp('fs', src)
      srcgt = 'ps';
    else
      srcgt = src;
    end
    
    xs = scenarios{ii,4};

    fS = aliasing_frequency(method, xs, src, 'position', xvec, conf);
    fS = reshape(fS, size(xS));
  
    if modus       
      suffix = sprintf(...
        'virtual %s, %s array (L = %d)\n%s', src, ...
        conf.secondary_sources.geometry, conf.secondary_sources.number, method);      
      
      % reference amplitude
      g = sound_field_mono(xref(1),xref(2),xref(3),[xs(1:3),0,1,0,1],srcgt,1,...
        f,conf);
      
      % compute synthesised sound field and aliasing error
      conf.resolution = 150;
      gtconf = conf;
      gtconf.secondary_sources.number = N0gt;
      gtconf.localwfs_vss.number = Nvgt;
      switch method
        case 'WFS'
          P = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);
          Pgt = sound_field_mono_wfs(X,Y,Z,xs,src,f,gtconf);
        case 'HOA'
          P = sound_field_mono_nfchoa(X,Y,Z,xs,src,f,conf);
          Pgt = sound_field_mono_nfchoa(X,Y,Z,xs,src,f,gtconf);
          suffix = sprintf('%s (M = %d)', suffix, conf.nfchoa.order);
        case 'LWFS-SBL'
          P = sound_field_mono_localwfs_sbl(X,Y,Z,xs,src,f,conf);
          Pgt = sound_field_mono_localwfs_sbl(X,Y,Z,xs,src,f,gtconf);
          suffix = sprintf('%s (M = %d, Npw = %d)', suffix, ...
            conf.localwfs_sbl.order, conf.localwfs_sbl.Npw);
        case 'LWFS-VSS'
          P = sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,conf);
          Pgt = sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,gtconf);

      end
      
      x0 = secondary_source_positions(conf);
      
      % plot aliasing frequency as function of x
      figure;
      subplot(2,2,1);
      imagesc(X,Y,fS)
      hold on
      contour(xS,yS,fS,[f,f],'k');
      hold off
      set(gca, 'YDir', 'normal', 'CLim', [f-1000 f+1000]);
      title(sprintf('aliasing frequency\n %s', suffix));
      
      % plot synthesised sound field
      subplot(2,2,2);
      imagesc(X,Y,real(P)./abs(g));
      hold on
      contour(xS,yS,fS,[f,f],'k');
      hold off
      set(gca, 'YDir', 'normal', 'CLim', [-1 1]);
      title(sprintf('synthesised sound field (f=%2.2f Hz)\n %s', f, suffix));
      
      % plot ground truth sound field
      subplot(2,2,3);
      imagesc(X,Y,real(Pgt)./abs(g));
      set(gca, 'YDir', 'normal', 'CLim', [-1 1]);
      title(sprintf('ground truth sound field (f=%2.2f Hz)\n %s', f, suffix));
      
      % plot aliasing error
      subplot(2,2,4);
      imagesc(X,Y,db((P-Pgt)./Pgt));
      hold on
      contour(xS,yS,fS,[f,f],'k');
      hold off
      set(gca, 'YDir', 'normal', 'CLim', [-40 10]);
      title(sprintf('aliasing error (f=%2.2f Hz)\n %s', f, suffix));
      
      for idx = 1:4
        subplot(2,2,idx);
        xlabel('x / m');
        ylabel('y / m');
        xlim(X);
        ylim(Y);
        hold on
        draw_loudspeakers(x0,[1 1 0], conf);
        hold off
        set_colorbar(conf);
        colorbar;
      end      
    end
end


status = true;
