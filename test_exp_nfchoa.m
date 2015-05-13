%% initialize
% close all;
clear variables;
% SFS Toolbox
SFS_start;

%% Parameters
conf = SFS_config_example;

conf.dimension = '2.5D';
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 60;
conf.secondary_sources.size = 3;
conf.secondary_sources.center = [0, 0, 0];

conf.plot.useplot = false;

conf.scattering.Nse = 23;
conf.scattering.Nce = 23;

conf.showprogress = true;
conf.resolution = 400;

ns = [0, -1, 0];  % propagation direction of plane wave
xs = [0,  2, 0];  % position of point source
f = 300;
xrange = [-2 2];
yrange = [-2 2];
zrange = 0;

xq = conf.secondary_sources.center;
xt = [ 0.5, 0.5, 0];
conf.xref = xt;
     
display(conf.scattering)

%% generic NFCHOA
% regular spherical expansion of plane wave and point source at xq
A1sph_original = sphexp_mono_pw(ns,f,xq,conf);
A2sph_original = sphexp_mono_ps(xs,'R', f,xq,conf);
% regular spherical expansion of plane wave and point source at xq+xt
A1sph = sphexp_mono_pw(ns,f,xq+xt,conf);
A2sph = sphexp_mono_ps(xs,'R', f,xq+xt,conf);
% regular-to-regular spherical reexpansion (translatory shift)
% [RRsph, RRsphm] = sphexp_mono_translation(-xt, 'RR', f, conf);
% shift spherical expansion back to xq
% A1sph_shift = RRsph*sphexp_bandlimit(A1sph,10);
% A2sph_shift = RRsph*sphexp_bandlimit(A2sph,10);
% loudspeakers
x0 = secondary_source_positions(conf);
% compute driving functions
D1sphhoa = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A1sph_original, f, conf);
D2sphhoa = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A2sph_original, f, conf);
% compute driving functions for shifted expansions
% D1sphhoa_shift = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A1sph_shift, f, conf);
% D2sphhoa_shift = driving_function_mono_nfchoa_sphexp(x0(:,1:3), A2sph_shift, f, conf);
% compute fields
P1sphhoa = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sphhoa,f,conf);
P2sphhoa = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sphhoa,f,conf);
% compute fields
% P1sphhoa_shift = sound_field_mono(xrange,yrange,zrange,x0,'ps',D1sphhoa_shift,f,conf);
% P2sphhoa_shift = sound_field_mono(xrange,yrange,zrange,x0,'ps',D2sphhoa_shift,f,conf);
% plot
[~,~,~,x1,y1,z1] = xyz_grid(xrange,yrange,zrange,conf);

plot_sound_field(P1sphhoa, x1,y1,z1, x0, conf);
title('NFCHOA: plane wave');
plot_scatterer(xq, norm(xt));
plot_sound_field(P2sphhoa ,x1,y1,z1, x0, conf);
title('NFCHOA: point source');
plot_scatterer(xq, norm(xt));
% plot_sound_field(P1sphhoa_shift/5, x1,y1,z1, x0, conf);
% title('NFCHOA: plane wave (shifted reexpansion)');
% plot_sound_field(P2sphhoa_shift/5, x1,y1,z1, x0, conf);
% title('NFCHOA: point source (shifted reexpansion)');

x1 = cosd(0:1:359).'*((0:0.1:1.4)*norm(xt)) + xq(1);
y1 = sind(0:1:359).'*((0:0.1:1.4)*norm(xt)) + xq(2);
z1 = zeros(size(x1)) + xq(3);


conf.showprogress = 0;
for rdx=1:size(x1,2)
  for idx=1:size(x1,1)  
    P1ring(idx,rdx) = sound_field_mono(x1(idx, rdx), y1(idx, rdx), z1(idx, rdx),x0,'ps',D1sphhoa,f,conf);
    P2ring(idx,rdx) = sound_field_mono(x1(idx, rdx), y1(idx, rdx), z1(idx, rdx),x0,'ps',D2sphhoa,f,conf);
    
    P1ring_gt(idx,rdx) = sound_field_mono(x1(idx, rdx), y1(idx, rdx), z1(idx, rdx),[ns, 1, 0, 0, 1], 'pw', 1, f,conf);
    P2ring_gt(idx,rdx) = sound_field_mono(x1(idx, rdx), y1(idx, rdx), z1(idx, rdx),[xs, 1, 0, 0, 1], 'ps', 1, f,conf);
  end
end
conf.showprogress = 1;

figure;
subplot(2,2,1);
imagesc((0:0.1:1.4), 0:1:359, db(P1ring./P1ring_gt));
colorbar;
colormap jet;
title('AMPLITUDE RATIO: plane wave');
subplot(2,2,2);
imagesc((0:0.1:1.4), 0:1:359, angle(P1ring./P1ring_gt));
colorbar;
colormap jet;
title('PHASE DIFFERENCE: plane wave');
subplot(2,2,3);
imagesc((0:0.1:1.4), 0:1:359, db(P2ring./P2ring_gt));
colorbar;
colormap jet;
title('AMPLITUDE RATIO: point source');
subplot(2,2,4);
imagesc((0:0.1:1.4), 0:1:359, angle(P2ring./P2ring_gt));
colorbar;
colormap jet;
title('PHASE DIFFERENCE: point source');
