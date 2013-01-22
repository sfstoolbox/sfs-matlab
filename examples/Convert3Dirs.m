%% Covert 3D Dataset
clear all
clc
close all

irs = load('az0.mat');

irs.apparent_azimuth = correct_azimuth(rad(irs.source_azimuth'));
irs.apparent_elevation = correct_elevation(rad(irs.source_elevation'));
irs = rmfield(irs,'source_azimuth');
irs = rmfield(irs,'source_elevation');
irs.ears = 'ear';
irs.head_position = [0;0;0];
irs.head_reference = [0;1;0];

irs.source_reference = [0;0;0];
irs.head_elevation = NaN;
irs.torso_azimuth = NaN;
irs.torso_elevation = NaN;
irs.distance = irs.distance';

irs = rmfield(irs,'microphone');
irs = rmfield(irs,'itd');


[x,y,z] = sph2cart(irs.apparent_azimuth,irs.apparent_elevation,irs.distance);
irs.source_position = [x;y;z];

save('FABIAN3D1,7.mat');