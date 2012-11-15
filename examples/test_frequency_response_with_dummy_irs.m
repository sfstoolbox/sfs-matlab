%% Test file to prove frequency response of a desired LS with dummy impulse
%  response (dirac) using 3D WFS
clear all
close all
clc
%% properties
conf = SFS_config;
L = 3;
conf.array = 'spherical';
conf.resolution_theta = 1;
X = [0 0];
phi = -pi/2;
src = 'pw';
xs = [0 2.5 0];
conf.xref = [0 0 0];
conf.usehpre = 1;
conf.usehcomp = 0;
dx1 = [];
dx2 = [];
%% get dirac impulses as impulse responses
irs = dummy_irs3d();
%% frequency response of the first LS
% 1. loudspeaker spacing = 0.16m
conf.dx0 = 0.16;
conf.hprefhigh = aliasing_frequency(conf.dx0);
ir = ir_wfs_3d(X,phi,xs,src,L,irs,conf);
[a,p,f] = easyfft(ir(:,1),conf);
dx1 = conf.dx0;
%% frequency response of the first LS
% 2. loudspeaker spacing = 0.32m
conf.dx0 = 0.32;
conf.hprefhigh = aliasing_frequency(conf.dx0);
ir = ir_wfs_3d(X,phi,xs,src,L,irs,conf);
[b,p,f] = easyfft(ir(:,1),conf);
dx2 = conf.dx0;
%% plot the results
figure; 
plot1 = semilogx(f,db(a),'b');
hold on
plot2 = semilogx(f,db(b),'r');
grid on
title('frequency response of LS 1/2')
xlabel('f/Hz')
ylabel('amplitude/dB')
legend([plot1 plot2],...
    sprintf('1st LS at %.2f m',dx1),...
    sprintf('2nd LS at %.2f m',dx2))
