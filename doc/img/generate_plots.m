%% ===== Secondary sources ===============================================
conf = SFS_config;
conf.secondary_sources.size = 3;

% === linear ===
conf.secondary_sources.geometry = 'line';
conf.secondary_sources.number = 21;
x0 = secondary_source_positions(conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0,conf);
axis([-2 2 -2 1]);
pause(1)
print_png('secondary_sources_linear.png');

% === circular ===
conf.secondary_sources.geometry = 'circle';
conf.secondary_sources.number = 56;
x0 = secondary_source_positions(conf);
figure;
figsize(540,404,'px');
draw_loudspeakers(x0,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('secondary_sources_circle.png');

% === box shaped ===
conf.secondary_sources.geometry = 'box';
conf.secondary_sources.number = 84;
x0 = secondary_source_positions(conf);
figure;
figsize(540,404,'px');
draw_loudspeakers(x0,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('secondary_sources_box.png');

% === box shaped with smoothed edges ===
conf.secondary_sources.geometry = 'rounded-box';
conf.secondary_sources.number = 84;
conf.secondary_sources.corner_radius = 0.3;
x0 = secondary_source_positions(conf);
figure;
figsize(540,404,'px');
draw_loudspeakers(x0,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('secondary_sources_rounded-box.png');

% === spherical array ===
conf.secondary_sources.geometry = 'sphere'; % or 'spherical'
conf.secondary_sources.number = 225;
x0 = secondary_source_positions(conf);
figure;
figsize(540,404,'px');
draw_loudspeakers(x0,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('secondary_sources_sphere.png');

% === arbitrary shaped arrays ===
% create a stadium like shape by combining two half circles with two linear
% arrays
% first getting a full circle with 56 loudspeakers
conf.secondary_sources.geometry = 'circle';
conf.secondary_sources.number = 56;
conf.secondary_sources.x0 = [];
x0 = secondary_source_positions(conf);
% store the first half cricle and move it up
x01 = x0(2:28,:);
x01(:,2) = x01(:,2) + ones(size(x01,1),1)*0.5;
% store the second half circle and move it down
x03 = x0(30:56,:);
x03(:,2) = x03(:,2) - ones(size(x03,1),1)*0.5;
% create a linear array
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 7;
conf.secondary_sources.size = 1;
x0 = secondary_source_positions(conf);
% rotate it and move it left
R = rotation_matrix(pi/2);
x02 = [(R*x0(:,1:3)')' (R*x0(:,4:6)')'];
x02(:,1) = x02(:,1) - ones(size(x0,1),1)*1.5;
x02(:,7) = x0(:,7);
% rotate it the other way around and move it right
R = rotation_matrix(-pi/2);
x04 = [(R*x0(:,1:3)')' (R*x0(:,4:6)')'];
x04(:,1) = x04(:,1) + ones(size(x0,1),1)*1.5;
x04(:,7) = x0(:,7);
% combine everything
conf.secondary_sources.geometry = 'custom';
conf.secondary_sources.x0 = [x01; x02; x03; x04];
% if we gave the conf.secondary_sources.x0 to the secondary_source_positions
% function it will simply return the defined x0 matrix
x0 = secondary_source_positions(conf);
figure;
figsize(540,404,'px');
draw_loudspeakers(x0,conf);
axis([-2 2 -2.5 2.5]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('secondary_sources_arbitrary.png');
conf.plot.realloudspeakers = true;
figure;
figsize(540,404,'px');
draw_loudspeakers(x0,conf);
axis([-2 2 -2.5 2.5]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('secondary_sources_arbitrary_realloudspeakers.png');

%% ===== Monochromatic sound fields ======================================
% === WFS 3D ===
conf = SFS_config;
conf.dimension = '3D';
conf.secondary_sources.size = 3;
conf.secondary_sources.number = 225;
conf.secondary_sources.geometry = 'sphere';
% [P,x,y,z,x0,win] = sound_field_mono_wfs_25d(X,Y,Z,xs,src,fconf);
sound_field_mono_wfs([-2 2],[-2 2],0,[0 -1 0],'pw',800,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_wfs_3d_xy.png');
sound_field_mono_wfs([-2 2],0,[-2 2],[0 -1 0],'pw',800,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_wfs_3d_xz.png');
sound_field_mono_wfs(0,[-2 2],[-2 2],[0 -1 0],'pw',800,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_wfs_3d_yz.png');
conf.resolution = 100;
sound_field_mono_wfs([-2 2],[-2 2],[-2 2],[0 -1 0],'pw',800,conf);
print_png('sound_field_wfs_3d_xyz.png');

% === WFS 2.5D ===
% simulating 2.5D WFS with circular array and a point source
conf = SFS_config;
conf.dimension = '2.5D';
conf.plot.useplot = true;
conf.plot.normalisation = 'center';
% [P,x,y,z,x0] = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);
[P,~,~,~,x0] = sound_field_mono_wfs([-2 2],[-2 2],0,[0 2.5 0],'ps',800,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_wfs_25d.png');
% plotting WFS with all secondary sources
x0_all = secondary_source_positions(conf);
[~,idx] = secondary_source_selection(x0_all,[0 2.5 0],'ps');
x0_all(:,7) = zeros(1,size(x0_all,1));
x0_all(idx,7) = x0(:,7);
conf.plot.realloudspeakers = true;
plot_sound_field(P,[-2 2],[-2 2],0,x0_all,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_wfs_25d_with_all_sources.png');
% simulating 2.5D NFCHOA with circular array and a plane wave
conf = SFS_config;
conf.dimension = '2.5D';
% sound_field_mono_nfchoa(X,Y,Z,xs,src,f,conf);
sound_field_mono_nfchoa([-2 2],[-2 2],0,[0 -1 0],'pw',800,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_nfchoa_25d.png');

% === 2D local WFS with box shaped array and circular virtual array ===
X = [-1 1];
Y = [-1 1];
Z = 0;
xs = [1 -1 0];
src = 'pw';
f = 7000;
conf = SFS_config;
conf.resolution = 1000;
conf.dimension = '2D';
conf.secondary_sources.geometry = 'box';
conf.secondary_sources.number = 4*56;
conf.secondary_sources.size = 2;
conf.localwfs_vss.size = 0.4;
conf.localwfs_vss.center = [0 0 0];
conf.localwfs_vss.geometry = 'circular';
conf.localwfs_vss.number = 56;
sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,conf);
axis([-1.1 1.1 -1.1 1.1]);
print_png('sound_field_localwfs_2d.png');

% === stereo setup ===
conf = SFS_config;
conf.plot.normalisation = 'center';
x0 = [-1 2 0 0 -1 0 1;1 2 0 0 -1 0 1];
% [P,x,y,z] = sound_field_mono(X,Y,Z,x0,src,D,f,conf)
sound_field_mono([-2 2],[-1 3],0,x0,'ps',[1 1],800,conf)
print_png('sound_field_stereo.png');

%% ===== spatio-temporal snapshots of the sound field ====================
conf = SFS_config;
conf.dimension = '2.5D';
conf.plot.useplot = true;
% sound_field_imp_nfchoa(X,Y,Z,xs,src,t,conf)
[p,x,y,z,x0] = sound_field_imp_nfchoa([-2 2],[-2 2],0,[0 2 0],'ps',0.005,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_imp_nfchoa_25d.png');
conf.plot.usedb = true;
plot_sound_field(p,[-2 2],[-2 2],0,x0,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_imp_nfchoa_25d_dB.png');
conf.plot.useplot = false;
conf.t0 = 'source';
t_40cm = 0.4/conf.c; % time to travel 40 cm in s
t0 = 0.0005; % start time of focused source in s
[p_ps,~,~,~,x0_ps] = ...
    sound_field_imp_wfs([-2 2],[-2 2],0,[1.9 0 0],'ps',t0+t_40cm,conf);
[p_pw,~,~,~,x0_pw] = ...
    sound_field_imp_wfs([-2 2],[-2 2],0,[1 -2 0],'pw',t0-t_40cm,conf);
[p_fs,~,~,~,x0_fs] = ...
    sound_field_imp_wfs([-2 2],[-2 2],0,[0 -1 0 0 1 0],'fs',t0,conf);
plot_sound_field(p_ps+p_pw+p_fs,[-2 2],[-2 2],0,[x0_ps; x0_pw; x0_fs],conf)
hold;
scatter(0,0,'k','x');   % origin of plane wave
scatter(1.9,0,'k','o'); % point source
scatter(0,-1,'k','o');  % focused source
hold off;
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_imp_multiple_sources_dB.png');

%% ===== custom grids ====================================================
conf = SFS_config;
conf.dimension = '3D';
conf.secondary_sources.number = 225;
conf.secondary_sources.geometry = 'sphere';
conf.resolution = 100;
conf.plot.normalisation = 'center';
X = randi([-2000 2000],125000,1)/1000;
Y = randi([-2000 2000],125000,1)/1000;
Z = randi([-2000 2000],125000,1)/1000;
sound_field_mono_wfs(X,Y,Z,[0 -1 0],'pw',800,conf);
print_png('sound_field_wfs_3d_xyz_custom_grid.png');
conf.plot.usedb = true;
conf.dimension = '2.5D';
conf.secondary_sources.number = 64;
conf.secondary_sources.geometry = 'circle';
sound_field_imp_nfchoa(X,Y,0,[0 2 0],'ps',0.005,conf);
axis([-2 2 -2 2]);
set(gca, 'XTick', -2:1:2); set(gca, 'YTick', -2:1:2);
print_png('sound_field_imp_nfchoa_25d_dB_custom_grid.png');

%% ===== modal windows ===================================================
conf = SFS_config;
conf.dimension = '2.5D';
conf.secondary_sources.number = 16;
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.size = 3;
conf.resolution = 300;
conf.plot.usedb = true;
conf.t0 = 'source';
X = [-2,2];
Y = [-2,2];
Z = 0;
conf.modal_window = 'rect';  % default
sound_field_imp_nfchoa(X,Y,Z,[0 -1 0],'pw',0,conf);
print_png('sound_field_imp_nfchoa_25d_dB_rect.png');
conf.modal_window = 'max-rE';
sound_field_imp_nfchoa(X,Y,Z,[0 -1 0],'pw',0,conf);
print_png('sound_field_imp_nfchoa_25d_dB_max-rE.png');
conf.modal_window = 'kaiser';
conf.modal_window_parameter = 1.0;
sound_field_imp_nfchoa(X,Y,Z,[0 -1 0],'pw',0,conf);
print_png('sound_field_imp_nfchoa_25d_dB_kaiser.png');
conf.modal_window = 'tukey';
conf.modal_window_parameter = 0.5;
sound_field_imp_nfchoa(X,Y,Z,[0 -1 0],'pw',0,conf);
print_png('sound_field_imp_nfchoa_25d_dB_tukey.png');


%% ===== impulse response of a spatial audio system ======================
conf = SFS_config;
conf.t0 = 'source';
X = [0 0 0];
phi = 0;
xs = [2.5 0 0];
src = 'ps';
t = (1:1000)/conf.fs*1000;
hrtf = dummy_irs(conf);
[ir,~,delay] = ir_wfs(X,phi,xs,src,hrtf,conf);
figure;
figsize(540,404,'px');
plot(t,ir(1:1000,1),'-g');
hold on;
offset = round(delay*conf.fs);
plot(t,ir(1+offset:1000+offset,1),'-b');
hold off;
xlabel('time / ms');
ylabel('amplitude');
print_png('impulse_response_wfs_25d.png');
X = [0 0 0];
head_orientation = [0 0];
xs = [2.5 0 0];
src = 'ps';
conf = SFS_config;
conf.N = 1000;
conf.t0 = 'source';
time_response_wfs(X,xs,src,conf)
axis([0 25 -0.005 0.025]);
print_png('impulse_response_wfs_25d_imp.png');


%% ===== frequency response of a spatial audio system ====================
X = [0 0 0];
head_orientation = [pi/2 0];
xs = [0 2.5 0];
src = 'ps';
hrtf = dummy_irs(conf);
conf = SFS_config;
conf.ir.usehcomp = false;
conf.wfs.usehpre = false;
[ir1,x0] = ir_wfs(X,head_orientation,xs,src,hrtf,conf);
conf.wfs.usehpre = true;
conf.wfs.hprefhigh = aliasing_frequency(x0,conf);
ir2 = ir_wfs(X,head_orientation,xs,src,hrtf,conf);
[a1,p,f] = spectrum_from_signal(norm_signal(ir1(:,1)),conf);
a2 = spectrum_from_signal(norm_signal(ir2(:,1)),conf);
figure;
figsize(540,404,'px');
semilogx(f,20*log10(a1),'-b',f,20*log10(a2),'-r');
axis([10 20000 -80 -40]);
set(gca,'XTick',[10 100 250 1000 5000 20000]);
legend('w/o pre-filter','w pre-filter');
xlabel('frequency / Hz');
ylabel('magnitude / dB');
print_png('frequency_response_wfs_25d.png');
% alternative variant
X = [0 0 0];
xs = [0 2.5 0];
src = 'ps';
conf = SFS_config;
freq_response_wfs(X,xs,src,conf);
axis([10 20000 -40 0]);
print_png('frequency_response_wfs_25d_mono.png');

