% generate dir for storing images
if ~exist('img','dir');
    mkdir('img');
end


%% ===== Secondary sources ===============================================
conf = SFS_config_example;
conf.secondary_sources.size = 3;

% === linear ===
conf.secondary_sources.geometry = 'line';
conf.secondary_sources.number = 21;
x0 = secondary_source_positions(conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0);
axis([-2 2 -2 1]);
pause(1)
print_png('img/secondary_sources_linear.png');

% === circular ===
conf.secondary_sources.geometry = 'circle';
conf.secondary_sources.number = 56;
x0 = secondary_source_positions(conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0);
axis([-2 2 -2 2]);
print_png('img/secondary_sources_circle.png');

% === box shaped ===
conf.secondary_sources.geometry = 'box';
conf.secondary_sources.number = 84;
x0 = secondary_source_positions(conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0);
axis([-2 2 -2 2]);
print_png('img/secondary_sources_box.png');

% === spherical array ===
conf.secondary_sources.geometry = 'sphere'; % or 'spherical'
conf.secondary_sources.number = 225;
x0 = secondary_source_positions(conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0);
axis([-2 2 -2 2]);
print_png('img/secondary_sources_sphere.png');

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
conf.secondary_sources.x0 = [x01; x02; x03; x04];
% if we gave the conf.secondary_sources.x0 to the secondary_source_positions
% function it will simply return the defined x0 matrix
x0 = secondary_source_positions(conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0);
axis([-2 2 -2.5 2.5]);
print_png('img/secondary_sources_arbitrary.png');
conf.plot.realloudspeakers = true;
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0,conf);
axis([-2 2 -2.5 2.5]);
print_png('img/secondary_sources_arbitrary_realloudspeakers.png');

%% ===== Monochromatic sound fields ======================================
% === stereo setup ===
conf = SFS_config_example;
x0 = [-1 2 0 0 -1 0 1;1 2 0 0 -1 0 1];
% [P,x,y,z] = sound_field_mono(X,Y,Z,x0,src,D,f,conf)
sound_field_mono([-2 2],[-1 3],0,x0,'ps',[1 1],800,conf)
print_png('img/sound_field_stereo.png');

% === WFS 3D ===
conf = SFS_config_example;
conf.dimension = '3D';
conf.secondary_sources.size = 3;
conf.secondary_sources.number = 225;
conf.secondary_sources.geometry = 'sphere';
% [P,x,y,z,x0,win] = sound_field_mono_wfs_25d(X,Y,Z,xs,src,fconf);
sound_field_mono_wfs([-2 2],[-2 2],0,[0 -1 0],'pw',800,conf);
print_png('img/sound_field_wfs_3d_xy.png');
sound_field_mono_wfs([-2 2],0,[-2 2],[0 -1 0],'pw',800,conf);
print_png('img/sound_field_wfs_3d_xz.png');
sound_field_mono_wfs(0,[-2 2],[-2 2],[0 -1 0],'pw',800,conf);
print_png('img/sound_field_wfs_3d_yz.png');


% simulating 2.5D WFS with circular array and a point source
conf = SFS_config_example;
conf.dimension = '2.5D';
conf.plot.useplot = 1;
% [P,x,y,z,x0] = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);
[P,x,y,z,x0] = sound_field_mono_wfs([-2 2],[-2 2],0,[0 2.5 0],'ps',800,conf);
print_png('img/sound_field_wfs_25d.png');
% plotting WFS with all secondary sources
x0_all = secondary_source_positions(conf);
[~,idx] = secondary_source_selection(x0,[0 2.5 0],'ps');
x0_all(:,7) = zeros(1,size(x0_all,1));
x0_all(idx,7) = x0(:,7);
plot_sound_field(P,x,y,z,x0_all,conf);
print_png('img/sound_field_wfs_25d_with_all_sources.png');
% simulating 2.5D NFCHOA with circular array and a plane wave
conf = SFS_config_example;
conf.dimension = '2.5D';
% sound_field_mono_nfchoa(X,Y,Z,xs,src,f,conf);
sound_field_mono_nfchoa([-2 2],[-2 2],0,[0 -1 0],'pw',800,conf);
print_png('img/sound_field_nfchoa_25d.png');

% --- spatio-temporal snapshots of the sound field ---
conf = SFS_config_example;
conf.dimension = '2.5D';
conf.plot.useplot = true;
% sound_field_imp_nfchoa(X,Y,Z,xs,src,t,conf)
[p,x,y,z,x0] = sound_field_imp_nfchoa([-2 2],[-2 2],0,[0 2 0],'ps',200,conf);
print_png('img/sound_field_imp_nfchoa_25d.png');
conf.plot.usedb = true;
conf.plot.colormap = 'jet';
plot_sound_field(p,x,y,z,x0,conf);
print_png('img/sound_field_imp_nfchoa_25d_dB.png');


% --- impulse response of the system ---
conf = SFS_config_example;
conf.ir.usehcomp = 0;
conf.wfs.usehpre = 0;
irs = dummy_irs;
[ir1,x0] = ir_wfs([0 0 0],pi/2,[0 2.5 0],'ps',irs,conf);
conf.wfs.usehpre = 1;
conf.wfs.hprefhigh = aliasing_frequency(x0);
ir2 = ir_wfs([0 0 0],pi/2,[0 2.5 0],'ps',irs,conf);
[a1,p,f] = easyfft(ir1(:,1)./max(abs(ir1(:,1))));
a2 = easyfft(ir2(:,1)./max(abs(ir2(:,1))));
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
semilogx(f,20*log10(a1),'-b',f,20*log10(a2),'-r');
axis([10 20000 -80 -40]);
set(gca,'XTick',[10 100 250 1000 5000 20000]);
legend('w/o pre-filter','w pre-filter');
xlabel('frequency / Hz');
ylabel('magnitude / dB');
print_png('img/impulse_response_wfs_25d.png');
% alternative variant
conf = SFS_config_example;
[a,f] = freq_response_wfs([0 0 0],[0 2.5 0],'ps',conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
semilogx(f,20*log10(a),'-r');
axis([10 20000 -20 20]);
set(gca,'XTick',[10 100 250 1000 5000 20000]);
legend('w pre-filter');
xlabel('frequency / Hz');
ylabel('magnitude / dB');
print_png('img/impulse_response_wfs_25d_mono.png');

% --- gnuplot ---
if ~system('gnuplot -V')
    conf = SFS_config_example;
    conf.plot.usegnuplot = 1;
    conf.plot.file = 'img/sound_field_nfchoa_25d_gnuplot.png';
    sound_field_mono_nfchoa([-2 2],[-2 2],0,[0 -1 0],'pw',1000,conf);
end
