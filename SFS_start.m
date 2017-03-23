function SFS_start(printbanner)
%SFS_START Start the Sound Field Synthesis Toolbox
%
%   Usage: SFS_start(printbanner)
%
%   Input parameters:
%       printbanner - 0: print nothing (default)
%                     1: print version and web page link
%
%   SFS_START(verbosity) starts the Sound Field Synthesis Toolbox (SFS).
%   This function must be run first in order to add the path's to Matlab.
%
%   See also: SFS_config, SFS_version

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


%% ===== Checking of input  parameters ===================================
nargmin = 0;
nargmax = 1;
narginchk(nargmin,nargmax);
if nargin==0
    printbanner = false;
end


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
    addpath([basepath,'/SFS_binaural_synthesis']);
    addpath([basepath,'/SFS_general']);
    addpath([basepath,'/SFS_helper']);
    addpath([basepath,'/SFS_ir']);
    addpath([basepath,'/SFS_monochromatic']);
    addpath([basepath,'/SFS_monochromatic/driving_functions_mono']);
    addpath([basepath,'/SFS_plotting']);
    addpath([basepath,'/SFS_plotting/colormaps']);
    addpath([basepath,'/SFS_ssr']);
    addpath([basepath,'/SFS_time_domain']);
    addpath([basepath,'/SFS_time_domain/driving_functions_imp']);
    addpath([basepath,'/SFS_time_domain/circexp_imp']);
    addpath([basepath,'/SFS_time_domain/pwd_imp']);
    addpath([basepath,'/SFS_HRTF_extrapolation']);
    addpath([basepath,'/validation']);
    if isoctave
        addpath([basepath,'/SFS_octave']);
    end
else
    path(path,basepath);
    path(path,[basepath,'/SFS_analysis']);
    path(path,[basepath,'/SFS_binaural_synthesis']);
    path(path,[basepath,'/SFS_general']);
    path(path,[basepath,'/SFS_helper']);
    path(path,[basepath,'/SFS_ir']);
    path(path,[basepath,'/SFS_monochromatic']);
    path(path,[basepath,'/SFS_monochromatic/driving_functions_mono']);
    path(path,[basepath,'/SFS_plotting']);
    path(path,[basepath,'/SFS_plotting/colormaps']);
    path(path,[basepath,'/SFS_ssr']);
    path(path,[basepath,'/SFS_time_domain']);
    path(path,[basepath,'/SFS_time_domain/driving_functions_imp']);
    path(path,[basepath,'/SFS_time_domain/circexp_imp']);
    path(path,[basepath,'/SFS_time_domain/pwd_imp']);
    path([basepath,'/SFS_HRTF_extrapolation']);
    path(path,[basepath,'/validation']);
    if isoctave
        path(path,[basepath,'/SFS_octave']);
    end
end


%% ===== Banner ==========================================================
if(printbanner)
    if ~usejava('desktop') % Looks only nice in console
        banner = sprintf( ...
            ['       ▄▄\n', ...
             ' ▄█▀▀ █▄ ▄█▀▀  Sound Field Synthesis Toolbox %s\n', ...
             ' ▄▄█▀ █  ▄▄█▀  http://matlab.sfstoolbox.org\n\n'], ...
            SFS_version);
    else
        banner = sprintf( ...
            ['\n', ...
             ' Sound Field Synthesis Toolbox %s\n', ...
             ' http://matlab.sfstoolbox.org\n\n'], ...
            SFS_version);
    end
    fprintf(banner);
end

