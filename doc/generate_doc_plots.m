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
% arbitrary shaped arrays
% create a stadium like shape by combining two half circles with two linear
% arrays
% first getting a full circle with 56 loudspeakers
conf.dx0 = L*pi/56;
conf.array = 'circle';
x0 = secondary_source_positions(L,conf);
% store the first half cricle and move it up
x01 = x0(2:28,:);
x01(:,2) = x01(:,2) + ones(size(x01,1),1)*0.5;
% store the second half circle and move it down
x03 = x0(30:56,:);
x03(:,2) = x03(:,2) - ones(size(x03,1),1)*0.5;
% create a linear array
conf.array = 'linear';
x0 = secondary_source_positions(1+conf.dx0,conf);
% rotate it and move it left
R = rotation_matrix(pi/2);
x02 = [(R*x0(:,1:2)')' x0(:,3) (R*x0(:,4:5)')' x0(:,6)];
x02(:,1) = x02(:,1) - ones(size(x0,1),1)*1.5;
% rotate it the other way around and move it right
R = rotation_matrix(-pi/2);
x04 = [(R*x0(:,1:2)')' x0(:,3) (R*x0(:,4:5)')' x0(:,6)];
x04(:,1) = x04(:,1) + ones(size(x0,1),1)*1.5;
% combine everything
conf.x0 = [x01; x02; x03; x04];
% if we gave the conf.x0 to the secondary_source_positions function it will
% simply return the defined x0 matrix
x0 = secondary_source_positions(L,conf);
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
draw_loudspeakers(x0);
axis([-2 2 -2.5 2.5]);
print_png('img/secondary_sources_arbitrary.png');

% --- monochromatic sound fields ---
% simulating stereo setup
conf = SFS_config_example;
% [x,y,P] = wave_field_mono_point_source(X,Y,xs,f);
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

% --- spatio-temporal snapshots of the sound field ---
conf = SFS_config_example;
conf.useplot = 1;
% wave_field_imp_nfchoa_25d(X,Y,xs,src,t,L,conf)
wave_field_imp_nfchoa_25d([-2 2],[-2 2],[0 2],'ps',200,3,conf);
print_png('img/wave_field_imp_nfchoa_25d.png');


% --- impulse response of the system ---
conf = SFS_config_example;
conf.usehcomp = 0;
conf.usehpre = 0;
irs = dummy_irs;
ir = ir_wfs_25d([0 0],pi/2,[0 2.5],'ps',3,irs,conf);
[a,p,f] = easyfft(ir(:,1)./max(abs(ir(:,1))));
figure;
figsize(conf.plot.size(1),conf.plot.size(2),conf.plot.size_unit);
semilogx(f,20*log10(a));
axis([10 20000 -100 -60]);
set(gca,'XTick',[10 100 250 1000 5000 20000]);
print_png('img/impulse_response_wfs_25d.png');


% --- gnuplot ---
conf = SFS_config_example;
conf.useplot = 1;
conf.plot.usegnuplot = 1;
conf.plot.file = 'img/wave_field_nfchoa_25d_gnuplot.png';
wave_field_mono_nfchoa_25d([-2 2],[-2 2],[0 -1],'pw',1000,3,conf);

