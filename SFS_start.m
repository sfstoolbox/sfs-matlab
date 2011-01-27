function SFS_start()
%SFS_START Start the Sound Field Synthesis Toolbox
%   Usage: SFS_start;
%
%   SFS_START starts the Sound Field Synthesis Toolbox (SFS). 
%   This function must be run first in order to add the path's to Matlab.
%
%   See also: SFS_config, SFS_version

% AUTHOR: Hagen Wierstorf


%% ===== Configuration ===================================================
printbanner = false;


%% ===== Adding Path's ===================================================

% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
basepath=which('SFS_start');
% Kill the function name from the path.
basepath=basepath(1:end-12);

% Add the base path and the needed sub-directories
if exist('addpath')
    addpath(basepath);
    addpath([basepath,'/SFS_analysis']);
    addpath([basepath,'/SFS_binaural_resynthesis']);
    addpath([basepath,'/SFS_general']);
    addpath([basepath,'/SFS_ir']);
    addpath([basepath,'/SFS_monochromatic']);
    addpath([basepath,'/SFS_plotting']);
    addpath([basepath,'/scripts']);
    addpath([basepath,'/SFS_time_domain']);
else
    path(path,basepath);
    path(path,[basepath,'/SFS_analysis']);
    path(path,[basepath,'/SFS_binaural_resynthesis']);
    path(path,[basepath,'/SFS_general']);
    path(path,[basepath,'/SFS_ir']);
    path(path,[basepath,'/SFS_monochromatic']);
    path(path,[basepath,'/SFS_plotting']);
    path(path,[basepath,'/scripts']);
    path(path,[basepath,'/SFS_time_domain']);
end


%% ===== Banner ==========================================================
if(printbanner)
    printf('SFS %1.1f successfully initialized.\n',SFS_version);
end

