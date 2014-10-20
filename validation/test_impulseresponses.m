%% initialize
close all;
clear variables;


%% Parameters
conf = SFS_config_example;
conf.showprogress = true;
conf.resolution = 400;
conf.plot.useplot = true;
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = false;
conf.plot.usedb = true;
conf.tapwinlen = 1.0;

% config for virtual array
conf.localsfs.method = 'wfs';
conf.localsfs.usetapwin = true;
conf.localsfs.vss.size = 1.0;
conf.localsfs.vss.center = [0, 0.5, 0];
conf.localsfs.vss.geometry = 'linear';
conf.localsfs.vss.number = 56;
conf.localsfs.vss.sampling = 'equi';
conf.localsfs.vss.logratio = 1.0;
conf.localsfs.vss.consider_target_field = true;
conf.localsfs.vss.consider_secondary_sources = true;
%conf.localsfs.vss.tapwinlen = 0.3;
%conf.localsfs.vss.wfs = conf.wfs;

% config for real array
conf.dimension = '2.5D';
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 56;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 1.5, 0];
conf.driving_functions = 'default';
conf.xref = conf.localsfs.vss.center;

xs = [0.0, -1.0, 0];  % propagation direction of plane wave
src = 'pw';

X = [-1.5 1.5];
Y = [-1, 1.55];
Z = 0;

%% temporal impulse responses
conf.ir.usehcomp = false;
irs = dummy_irs;

s_lwfs = ir_localwfs(conf.xref,pi/2,xs,src,irs,conf);
s_wfs = ir_wfs(conf.xref,pi/2,xs,src,irs,conf);

[S_lwfs, ~, f_lwfs] = easyfft(s_lwfs(:,1)./max(abs(s_lwfs(:,1))), conf);
[S_wfs, ~, f_wfs] = easyfft(s_wfs(:,1)./max(abs(s_wfs(:,1))), conf);

%% spatio-temporal sound field
sound_field_imp_localwfs(X,Y,Z, xs, src, 400, conf);
sound_field_imp_wfs(X,Y,Z, xs, src, 190, conf);
