function status = test_imp_25d(modus)
%TEST_imp_25D tests behavior of 2.5D SFS techniques in time-domain
%
%   Usage: status = test_imp_25d(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false

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
conf.xref = [0,0,0];
conf.dimension = '2.5D';
conf.plot.useplot = false;
conf.plot.usedb = true;
conf.wfs.hpreflow = 20;
conf.wfs.hprefhigh = 20000;
conf.localwfs_vss.wfs = conf.wfs;
conf.localwfs_vss.geometry = 'circular';
conf.localwfs_vss.number = 512;
conf.localwfs_vss.center = conf.xref;
conf.localwfs_vss.size = 0.6;
conf.localwfs_vss.consider_secondary_sources = false;
conf.localwfs_vss.consider_target_field = false;
conf.localwfs_sbl.order = 23;
conf.delayline.resamplingfactor = 8;
conf.delayline.resampling = 'pm';
conf.delayline.filter = 'lagrange';
conf.delayline.filterorder = 9;
conf.t0 = 'source';

X = [-2.05 2.05];
Y = [-2.05 2.05];
Z = 0;

% test scenarios
scenarios = { ...
  'WFS', 'reference_point', 'linear', 'pw', [ 0.0 -1.0   0.0]
  'WFS', 'reference_point', 'linear', 'ps', [ 0.0  2.5   0.0]
  'WFS', 'reference_point', 'linear', 'fs', [ 0.0  0.75  0.0  0.0 -1.0  0.0]
  'WFS', 'reference_line' , 'linear', 'pw', [ 0.0 -1.0   0.0]
  'WFS', 'reference_line' , 'linear', 'ps', [ 0.0  2.5   0.0]
  'WFS', 'reference_line' , 'linear', 'fs', [ 0.0  0.75  0.0  0.0 -1.0  0.0]
  'HOA', 'default', 'circular', 'pw', [ 0.0 -1.0   0.0]
  'HOA', 'default', 'circular', 'ps', [ 0.0  2.5  0.0]
  'LWFS-VSS', 'default', 'circular', 'pw', [ 0.0 -1.0   0.0]
  'LWFS-VSS', 'default', 'circular', 'ps', [ 0.0  2.5  0.0]
  'LWFS-SBL', 'default', 'circular', 'pw', [ 0.0 -1.0   0.0]
  'LWFS-SBL', 'default', 'circular', 'ps', [ 0.0  2.5  0.0]
  };

%% ===== Main ============================================================

sofa = dummy_irs(512,conf);
xt = conf.xref + [0.1, 0, 0];

for ii=1:size(scenarios)
  
    src = scenarios{ii,4};  %
    xs = scenarios{ii,5};  % source position
    switch src
    case 'pw'
        tau = xs*conf.xref.'./conf.c;
    case {'ps', 'fs'}
        tau = norm(xs(1:3) - conf.xref)./conf.c;
    end
    
    % get listening area
    conf.secondary_sources.geometry = scenarios{ii,3};
    switch scenarios{ii,3}
    case 'linear'
        conf.secondary_sources.size = 4;
        conf.secondary_sources.number = 128;
        conf.usetapwin = true;
        conf.tapwinlen = 0.2;
        conf.secondary_sources.center = [0, 1.5, 0];
    case 'circular'
        conf.secondary_sources.size = 3;
        conf.secondary_sources.number = 128;
        conf.secondary_sources.center = [0, 0, 0];
    end
    x0 = secondary_source_positions(conf);
    
    % compute driving signals
    conf.driving_functions = scenarios{ii,2};
    switch scenarios{ii,1}
    case 'WFS'
        x0 = secondary_source_selection(x0,xs,src);
        x0 = secondary_source_tapering(x0,conf);
        [d,~,~,delay_offset] = driving_function_imp_wfs(x0,xs,src,conf);
    case 'HOA'
        [d,~,delay_offset] = driving_function_imp_nfchoa(x0,xs,src,conf);
    case 'LWFS-VSS'
        [d,x0,~,~,delay_offset] = ...
            driving_function_imp_localwfs_vss(x0,xs,src,conf);
    case 'LWFS-SBL'
        [d,delay_offset] = driving_function_imp_localwfs_sbl(x0,xs,src,conf);
    end
    
    % Spatio-temporal sound field
    p = sound_field_imp(X,Y,Z,x0,'ps',d,tau+delay_offset,conf);
    
    % spectrum of reproduced sound field at reference position
    ir_sfs = ir_generic(xt,0,x0,d,sofa,conf);
    [IR_sfs,~,f_sfs] = spectrum_from_signal(ir_sfs(:,1),conf);
    
    % spectrum of ground truth sound field at reference position
    if strcmp(src,'pw')
        ir_gt = ir_point_source(xt,0,-xs./norm(xs),sofa,conf);
        ir_gt = ir_gt*4*pi;
    else
        ir_gt = ir_point_source(xt,0,xs(1:3),sofa,conf);
    end
    [IR_gt,~,f_gt] = spectrum_from_signal(ir_gt(:,1),conf);
    
    if modus
        plot_sound_field(p,X,Y,Z,x0,conf);
        title(sprintf('%s %s %s',scenarios{ii,1},src, ...
                      conf.driving_functions),'Interpreter','none');
        figure;
        semilogx(f_sfs,db(IR_sfs),'r',f_gt,db(IR_gt),'b--');
        xlabel('Frequency / Hz');
        ylabel('Magnitude / dB');
        title(sprintf('%s %s %s',scenarios{ii,1},src,conf.driving_functions), ...
            'Interpreter','none');
        legend('reproduced','ground truth','Location','northwest');
    end  
end

status = true;
