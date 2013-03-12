%%FIND_ARC_OF_CIRCLES_FOR_DIFFERENT_ELEVATION_ANGLE finds arc of cicles for 
% different elevation angles of a 3D HRIR dataset
clear all;close all;clc

load('irs_pw.mat');


%% 
idx = [];
output_save_name = 'irs_pw_theta_%d';
output_folder_path = 'irs_splitted/';

for ii=2:length(irs.apparent_elevation)-1;
    
    idx = find(irs_pw.apparent_elevation == irs_pw.apparent_elevation(1,ii)); 
%     save_file_name = num2str(degree(irs.apparent_elevation(ii)),'%02d');
    save_file_name = sprintf(output_save_name,degree(irs_pw.apparent_elevation(ii)));
    ii = idx(end);
    
    irs_tmp.apparent_azimuth = irs_pw.apparent_azimuth(idx);
    irs_tmp.apparent_elevation = irs_pw.apparent_elevation(idx);
    irs_tmp.distance = 'inf';%irs_pw.distance(idx);
    irs_tmp.left = irs_pw.left(:,idx);
    irs_tmp.right = irs_pw.right(:,idx);
    save([output_folder_path,save_file_name,'.mat'],'-struct','irs_tmp');
end