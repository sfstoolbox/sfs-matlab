%TEST_EXTRAPOLATE_FARFIELD_HRTF_25D tests the far field range extrapolation
% of a given 2.5d HRTF dataset. It is necessary to have the dataset as irs-
% structure in the folder where this function is located.
% After the range extrapolation the ILD of the irs dataset and the
% extrapolated dataset is compared.

%% ===== Configuration ===================================================
clear all; close all; clc;
conf = SFS_config;
conf.array = 'circle';
irs = read_irs('QU_KEMAR_anechoic_3m.mat');
R = irs.distance;
L = 2*R;
nls = length(irs.apparent_azimuth);
% potential error sources
conf.usefracdelay = 0;
conf.usetapwin = 1;
conf.usehpre = 1;
%% ===== Computation =====================================================
% far field extrapolation
irs_pw = extrapolate_farfield_hrtfset_25d(irs,conf);

%calculated ILD
ild1 = interaural_level_difference(irs.left(:,:),irs.right(:,:));
ild2 = interaural_level_difference(irs_pw.left(:,:),irs_pw.right(:,:));
%% ===== Plot ============================================================
% ILD
figure
plot(ild1,'r')
hold on
plot(ild2,'b')
grid on
legend('without extrapolation','with extrapolation')
title('ILD')