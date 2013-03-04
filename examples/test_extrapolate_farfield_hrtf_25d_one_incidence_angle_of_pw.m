%TEST_EXTRAPOLATE_FARFIELD_HRTF_25D_ONE_INCIDENCE_ANGLE_OF_PW test the
% range extrapolation of ONE desired incidence angle of a plane wave(2.5d).
% The incidence angle can be choosen by 'position'.
% After calculating the frequency response of the extrapolated HRIR and the
% measured HRIR will be plotted.
% The used HRIR dataset has to be stored in the same folder where this
% function is located.

%% ===== Configuration ===================================================
clear all; close all; clc;
conf = SFS_config;
conf.array = 'circle';
irs = read_irs('QU_KEMAR_anechoic_3m.mat');
R = irs.distance;
L = 2*R;
nls = length(irs.apparent_azimuth);

% potential error sources
conf.usefracdelay = 1;
conf.usetapwin = 1;
conf.usehpre = 1;

% Desired HRTF position at azimuth
position = -90; % in degree

idx = findrows(...
    [irs.apparent_azimuth' irs.apparent_elevation'],...
    [rad(position),rad(0)]);

%% check the direction first
ir = get_ir(irs,rad(position),rad(0),R,conf.xref);
% t = 0:1/44100:10;
% test = auralize_ir(ir,sin(2*pi*440*t),1,conf);
% wavwrite(test,conf.fs,32,'test');
% x = sin(2*pi*440*0:1/44100:1/44100*10);
% outsig = auralize_ir(ir,x,1,conf);
% wavwrite(outsig,conf.fs,32,'HRTF');
fs = conf.fs;                   % sampling frequency
%% ===== Variables ======================================================

Acorr = -1.7;                     % DAGA 2011 R=0.5m -> pw
Af = Acorr*sin(irs.apparent_azimuth(idx));

%% HRTFs:       
%               phi = 0 degree --> xs = [0 -1 0] direction of pw
%             $$$$$$$$$$$$$$$$$
%           $                   $
%         $                        $ 
%      $                             $                  y
%    $                                 $                ^
%   $                                    $              |
%  $                                      $             |
%  $                                       $            |_ _ _ _>  x
% $                                         $                                 
% $                    o                    $ -90 degree
% $            listener position            $
% $                 [0 0 0]                 $
% $                                        $
%  $                                      $
%   $                                    $
%    $                                  $
%      $                               $
%       $                           $
%          $                     $
%              $$$$$$$$$$$$$$$$
%                +/- 180 degree
%          
%          
% get virtual loudspeaker positions from HRTF dataset
conf.array = 'circle';
conf.dx0 = 2*pi*R/nls;
conf.hprefhigh = aliasing_frequency(conf.dx0,conf);
conf.xref = [0 0 0];
conf.x0 = zeros(nls,6);
conf.x0(:,1:2) = [R*cos(irs.apparent_azimuth) ; R*sin(irs.apparent_azimuth)]';
conf.x0(:,4:6) = direction_vector(conf.x0(:,1:3),repmat(conf.xref,nls,1));


%% ===== Computation =====================================================
% get virtual secondary source positions
x0_help = secondary_source_positions(L,conf);

% Initialize new irs set
irs_pw = irs;
irs_pw.description = 'Extrapolated HRTF set containing plane waves';
irs_pw.left = zeros(length(irs_pw.left),1);
irs_pw.right = zeros(length(irs_pw.right),1);
irs_pw.distance = 'Inf';

% Generate a irs set for all given angles


    % direction of plane wave
%     xs = -[cos(rad(position)) sin(rad(position)) 0];
    xs = [-1 0 0];

    % calculate active virtual speakers
    x0 = secondary_source_selection(x0_help,xs,'pw');

    
    % generate tapering window
    win = tapering_window(x0,conf);

    % sum up contributions from individual virtual speakers
    for l=1:size(x0,1)
        % Driving function to get weighting and delaying
        [a,delay] = driving_function_imp_wfs_25d(x0(l,:),xs,'pw',conf);
        dt = delay*fs + round(R/conf.c*fs);
        w=a*win(l);
        % truncate IR length
        irl = fix_ir_length(irs.left(:,l),length(irs.left(:,l)),0);
        irr = fix_ir_length(irs.right(:,l),length(irs.right(:,l)),0);
        % delay and weight HRTFs
        irs_pw.left(:,1) = irs_pw.left(:,1) + delayline(irl',dt,w,conf)';
        irs_pw.right(:,1) = irs_pw.right(:,1) + delayline(irr',dt,w,conf)';
    end

    irs_pw.left(:,1) = irs_pw.left(:,1)*10^(Af/20);
    irs_pw.right(:,1) = irs_pw.right(:,1)*10^(-Af/20);



%% ===== Pre-equalization ===============================================
irs_pw.left = wfs_preequalization(irs_pw.left(:,1),conf);
irs_pw.right = wfs_preequalization(irs_pw.right(:,1),conf);


%% plot the stuff

% plot wavefield with active secondary sources in frequency domain
conf.zreferenceaxis = 'y';
conf.useplot = 1;
wave_field_mono_wfs_25d([-4 4],[-4 4],[0 0],xs,'pw',1000,L,conf);

%calculated ILD
ild1 = interaural_level_difference(ir(:,:),ir(:,:));
ild2 = interaural_level_difference(irs_pw.left(:,:),irs_pw.right(:,:));

% plot ILD
figure
plot(ild1,'rx','MarkerSize',10)
hold on
plot(ild2,'bx','MarkerSize',10)
grid on
legend('without extrapolation','with extrapolation')
title('ILD')

% FFT of HRIRs and extrapolated HRIRs
[amplitude1,phase1,f] = easyfft(ir(:,1),conf);
[amplitude2,phase2,f] = easyfft(irs_pw.left(:,1),conf);
[amplitude3,phase3,f] = easyfft(ir(:,2),conf);
[amplitude4,phase4,f] = easyfft(irs_pw.right(:,1),conf);

% plot frequency response left ear
figure
plot(f,db(amplitude1),'r')
hold on
plot(f,db(amplitude2),'b')
grid on
legend('without extrapolation','with extrapolation')
title('left ear')

% plot frequency response right ear
figure
plot(f,db(amplitude3),'r')
hold on
plot(f,db(amplitude4),'b')
grid on
legend('without extrapolation','with extrapolation')
title('right ear')

