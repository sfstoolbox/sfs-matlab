% generate dir for storing images
if ~exist('img','dir');
    mkdir('img');
end

% --- secondary sources ---
conf = SFS_config_example;
L = 3;
% linear
conf.array = 'linear';
x0 = secondary_source_positions(L,conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0);
axis([-2 2 -2 1]);
print_png('img/secondary_sources_linear.png');
% circular
conf.array = 'circle';
x0 = secondary_source_positions(L,conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0);
axis([-2 2 -2 2]);
print_png('img/secondary_sources_circle.png');
% box shaped
conf.array = 'box';
x0 = secondary_source_positions(L,conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0);
axis([-2 2 -2 2]);
print_png('img/secondary_sources_box.png');

% --- monochromatic sound fields ---
% simulating stereo setup
% [x,y,P] = wave_field_mono_point_source(X,Y,xs,f);
conf = SFS_config_example;
[x,y,P1] = wave_field_mono_point_source([-2 2],[-1 3],[-1 2],1000);
[x,y,P2] = wave_field_mono_point_source([-2 2],[-1 3],[1 2],1000);
plot_wavefield(x,y,real(P1+P2),[-1 2 0 0 -1 0;1 2 0 0 -1 0],conf);
print_png('img/wave_field_stereo.png');
% simulating 2.5D WFS with circular array and a point source
conf = SFS_config_example;
conf.useplot = 1;
% [x,y,P,x0,win] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
[x,y,P,x0,win] = wave_field_mono_wfs_25d([-2 2],[-2 2],[0 2.5],'ps',1000,3,conf);
print_png('img/wave_field_wfs_25d.png');
% plotting WFS with all secondary sources
x0 = secondary_source_positions(L,conf);
[~,idx] = secondary_source_selection(x0,[0 2.5],'ps');
win2 = zeros(1,size(x0,1));
win2(idx) = win;
plot_wavefield(x,y,P,x0,win2,conf);
print_png('img/wave_field_wfs_25d_with_all_sources.png');
% simulating 2.5D NFCHOA with circular array and a plane wave
conf = SFS_config_example;
conf.useplot = 1;
% wave_field_mono_nfchoa_25d(X,Y,xs,src,f,L,conf);
wave_field_mono_nfchoa_25d([-2 2],[-2 2],[0 -1],'pw',1000,3,conf);
print_png('img/wave_field_nfchoa_25d.png');
