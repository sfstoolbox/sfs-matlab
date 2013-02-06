%% Test auralization 3D
clear all; close all; clc;
conf = SFS_config;
conf.array = 'spherical';
conf.usehpre = 1;
conf.fs = 44100;
irs = read_irs('FABIAN_3D_anechoic_~1.6.mat');

%% Desired HRTF position at azimuth

position = 90; % in degree

idx = findrows(...
    [irs.apparent_azimuth' irs.apparent_elevation'],...
    [rad(position),rad(0)]);
R = irs.distance(1,idx);

%% check the direction first
ir = get_ir(irs,rad(position),rad(0),R,[0 0 0]);
x = wavread('Sinus440Hz');
outsig = auralize_ir(ir,x,1,conf);
wavwrite(outsig,conf.fs,32,'HRTF');

%%
fs = conf.fs;
L = 2*R;
conf.dx0 = R*pi/180*2;
conf.hprefhigh = aliasing_frequency(conf.dx0,conf);
conf.xref = [0 0 0];
x0(:,1:3) = irs.source_position.';
x0(:,4:6) = direction_vector(x0(:,1:3),repmat(conf.xref,length(irs.left),1));
x0(:,8) = ones(length(irs.apparent_elevation),1);%sin(irs.apparent_elevation);
x0(:,7) = ones(length(irs.apparent_elevation),1);

% Initialize new irs set
irs_pw = irs;
irs_pw.description = 'Extrapolated HRTF set containing plane waves';
irs_pw.left = zeros(size(irs_pw.left));
irs_pw.right = zeros(size(irs_pw.right));
irs_pw.distance = 'Inf';

% direction of plane wave
xs = -1*[cos(rad(position)) sin(rad(position)) 0];

% calculate active virtual speakers
x0 = secondary_source_selection(x0,xs,'pw');
    
% sum up contributions from individual virtual speakers
for l=1:size(x0,1)
    % Driving function to get weighting and delaying
    [w,delay] = driving_function_imp_wfs_3d(x0(l,:),xs,'pw',conf);
    dt = delay*fs + round(norm(x0(l,1:3))/conf.c*fs);  
    % truncate IR length
    irl = fix_ir_length(irs.left(:,l),length(irs.left(:,l)),0);
    irr = fix_ir_length(irs.right(:,l),length(irs.right(:,l)),0);
    % delay and weight HRTFs
    irs_pw.left(:,1) = irs_pw.left(:,1) + delayline(irl',dt,w,conf)';
    irs_pw.right(:,1) = irs_pw.right(:,1) + delayline(irr',dt,w,conf)';
end


%% ===== Pre-equalization ===============================================
irs_pw.left = wfs_preequalization3d(irs_pw.left,conf);
irs_pw.right = wfs_preequalization3d(irs_pw.right,conf);

outsig = auralize_ir([irs_pw.left(:,1) irs_pw.right(:,1)],x,1,conf);
wavwrite(outsig,conf.fs,32,'HRTF2');

%% plot the stuff
ild1 = interaural_level_difference(ir(:,1),ir(:,2));
ild2 = interaural_level_difference(irs_pw.left(:,1),irs_pw.right(:,1));

plot(ild1,'rx','MarkerSize',10)
hold on
plot(ild2,'bx','MarkerSize',10)
grid on

[amplitude1,phase1,f] = easyfft(ir(:,1),conf);
[amplitude2,phase2,f] = easyfft(irs_pw.left(:,1),conf);

figure
plot(f,db(amplitude1),'r','MarkerSize',10)
hold on
plot(f,db(amplitude2),'b','MarkerSize',10)
grid on

%%  Show wavefield in the Frequency domain
%% properties of SFS_conf
% 
% conf = SFS_config;
% conf.zreferenceaxis = 'y'; % set a plot reference axis
% conf.useplot = 1; % plot the results = 1 ; otherwise = 0
% conf.array = 'spherical'; % array type
% conf.usetapwin = 0; % do not use tapering window, it's isn't needed in the 3D case
% conf.plot.loudspeakers = 0; % do not plot loudspeakers in the 3D case, because it's a mess ;)
% conf.usehpre = 1; % use preequalization filter
% % conf.number_of_points_on_sphere = 81^2; % number of points on the sphere, if spherical array is choosen
% conf.debug = 1; % debug=1 allows to plot results of different evaluation steps
% conf.frame = 0; % set a frame to show the wavefield in the time domain
% conf.xref = [0 0 0]; % ps: 'listener position' ; pw:  place where the wavefield is scaled to one
% conf.grid = 'HRTFgrid';
% %% properties of desired wavefield
% % xs = -[1 0 0]; % position of point source / inicidence angle of plane wave
% r = 0; % radius of the sphere
% L = 2*r; % diameter of the sphere
% src = 'pw'; % select source type pw/ps/fs
% f = 1000; % evaluation frequency
% %% plot properties
% scale_axis_1 = [0 0];
% scale_axis_2 = [-2 2];
% scale_axis_3 = [0.5 0.5];
% 
% %% calculate frequency domain wavefield
% wave_field_mono_wfs_3d(scale_axis_2,scale_axis_2,scale_axis_1,xs,src,f,L,conf);
% wave_field_imp_wfs_3d(scale_axis_2,scale_axis_2,scale_axis_1,xs,src,L,conf);
% grid on
