function SFS_stop()
%SFS_STOP Stop the Sound Field Synthesis Toolbox
%
%   Usage: SFS_stop;
%
%   SFS_STOP stops the Sound Field Synthesis Toolbox (SFS). 
%   This function removes the path's of the toolbox.
%
%   see also: SFS_start, SFS_config, SFS_version

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Configuration ===================================================
printbanner = false;


%% ===== Adding Path's ===================================================

% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
basepath=which('SFS_stop');
% Kill the function name from the path.
basepath=basepath(1:end-11);

% Add the base path and the needed sub-directories
if exist('rmpath')
    if isoctave
        rmpath([basepath,'/SFS_octave']);
    end
    
    rmpath([basepath,'/SFS_analysis']);
    rmpath([basepath,'/SFS_binaural_synthesis']);
    rmpath([basepath,'/SFS_general']);
    rmpath([basepath,'/SFS_helper']);
    rmpath([basepath,'/SFS_ir']);
    rmpath([basepath,'/SFS_monochromatic']);
    rmpath([basepath,'/SFS_monochromatic/driving_functions_mono']);
    rmpath([basepath,'/SFS_plotting']);
    rmpath([basepath,'/SFS_ssr']);
    rmpath([basepath,'/SFS_time_domain']);
    rmpath([basepath,'/SFS_time_domain/driving_functions_imp']);
    rmpath([basepath,'/SFS_HRTF_extrapolation']);
    rmpath([basepath,'/validation']);
    %addpath(basepath);
end


%% ===== Banner ==========================================================
if(printbanner)
    printf('SFS %1.1f successfully stopped.\n',SFS_version);
end

