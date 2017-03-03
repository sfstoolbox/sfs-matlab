function SFS_stop()
%SFS_STOP Stop the Sound Field Synthesis Toolbox
%
%   Usage: SFS_stop;
%
%   SFS_STOP stops the Sound Field Synthesis Toolbox (SFS). 
%   This function removes the path's of the toolbox.
%
%   See also: SFS_start, SFS_config, SFS_version

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
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
    rmpath([basepath,'/SFS_plotting/colormaps']);
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

