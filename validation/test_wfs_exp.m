function boolean = test_wfs_exp
%TEST_WFS_EXP tests the correctness of the driving functions

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

% TODO: add mode to save data as reference data


%% ===== Checking of input  parameters ===================================
nargmin = 0;
nargmax = 0;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
conf = SFS_config;
conf.secondary_sources.size = 3;
f = 1000;
conf.plot.useplot = false;
conf.plot.usenormalisation = true;
conf.plot.normalisation = 'center';
conf.driving_functions = 'default';
conf.usetapwin = false;
Nce = 15;

% test scenarios
scenarios = { ...
    'WFS', '2D',   'circular', 'pw', [ 0.5  0.5  0.0], 'circ'; ...
    'WFS', '2D',   'circular', 'ls', [ 0.0  2.5  0.0], 'circ'; ...
    'WFS', '2.5D', 'circular', 'pw', [ 0.5  0.5  0.0], 'circ'; ...
    'WFS', '2.5D', 'circular', 'pw', [ 0.5  0.5  0.0], 'sph'; ...
    'WFS', '2.5D', 'circular', 'ls', [ 0.0  2.5  0.0], 'circ'; ...
    'WFS', '2.5D', 'circular', 'ls', [ 0.0  2.5  0.0], 'sph'; ...
    'WFS', '2.5D', 'circular', 'ps', [ 0.0  2.5  0.0], 'sph'; ...
    'WFS', '3D', 'sphere', 'pw', [ 0.5  0.5  0.0], 'sph'; ...
    'WFS', '3D', 'sphere', 'ls', [ 0.0  2.5  0.0], 'sph'; ...
    'WFS', '3D', 'sphere', 'ps', [ 0.0  2.5  0.0], 'sph'; ...
};

% Start testing
for ii=1:size(scenarios)

    % get current dimension
    conf.dimension = scenarios{ii,2};

    % get source type and position
    src = scenarios{ii,4};
    xs = scenarios{ii,5};
    
    % get secondary source type
    switch conf.dimension
      case '2D'
        secsrc = 'ls';
      case {'2.5D', '3D'}
        secsrc = 'ps';
    end
    
    % get secondary source distribution
    conf.secondary_sources.geometry = scenarios{ii,3};
    switch conf.secondary_sources.geometry
        case 'linear'
            X = [-2 2];
            Y = [-3 0.15];
            Z = 0;
            conf.xref = [0 -1.5 0];
            conf.secondary_sources.number = 20;
        case 'circular'
            X = [-2 2];
            Y = [-2 2];
            Z = 0;
            conf.xref = [0 0 0];
            conf.secondary_sources.number = 56;
        case 'box'
            X = [-2 2];
            Y = [-2 2];
            Z = 0;
            conf.xref = [0 0 0];
            conf.secondary_sources.number = 80;
        case 'sphere'
            X = [-2 2];
            Y = [-2 2];
            Z = 0;
            conf.xref = [0 0 0];
            conf.secondary_sources.number = 900;
    end
    x0 = secondary_source_positions(conf);
    
    % get expansion coefficients
    xq = conf.xref;
    fname = [scenarios{ii,6}, 'exp_mono_', src];
    expfunc = str2func(fname);
    switch src
      case 'pw'
        A = expfunc(xs, Nce,f,xq,conf);
      case 'ls' 
        A = expfunc(xs,'R', Nce,f,xq,conf);
      case 'ps'
        A = expfunc(xs,'R', Nce,f,xq,conf);
    end    
    
    % get method
    method = scenarios{ii,1};
    
    % ===== WFS ==========================================================
    if strcmp('WFS',method)
      x0 = secondary_source_selection(x0,xs,src);
      x0 = secondary_source_tapering(x0,conf);      
      
      Dfunc = str2func(['driving_function_mono_wfs_', scenarios{ii,6}, 'exp']); 
        
      D = Dfunc(x0(:,1:3),x0(:,4:6),A,'R',f,xq,conf);
      P = sound_field_mono(X,Y,Z,x0,secsrc,D,f,conf);
      plot_sound_field(P, X, Y, Z, x0, conf);
    end
    
    title(sprintf('%s %s of %s. exp. of %s', conf.dimension, method, ...
      scenarios{ii,6}, src), 'Interpreter', 'None');    
end
