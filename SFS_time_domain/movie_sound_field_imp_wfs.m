function movie_sound_field_imp_wfs(X,Y,Z,xs,src,outfile,conf)
%MOVIE_SOUND_FIELD_IMP_WFS generates movie a WFS sound field
%
%   Usage: movie_sound_field_imp_wfs(X,Y,Z,xs,src,outfile,conf)
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax] or nD-array 
%       Y           - y-axis / m; single value or [ymin,ymax] or nD-array
%       Z           - z-axis / m; single value or [zmin,zmax] or nD-array
%       xs          - position of point source / m
%       src         - sourcetype of the virtual source:
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       outfile     - name for the movie file
%       conf        - configuration struct (see SFS_config)
%
%   MOVIE_SOUND_FIELD_IMP_WFS(X,Y,Z,xs,src,outfile,conf) generates a movie of
%   simulations of a sound field of the given source positioned at xs
%   using a WFS driving function in the temporal domain with different phase.
%
%   See also: sound_field_imp_wfs, plot_sound_field, generate_movie

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


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 7;
narginchk(nargmin,nargmax);
isargxs(xs);
isargchar(src,outfile);
isargstruct(conf);


%% ===== Configuration ==================================================
% Plotting
useplot = conf.plot.useplot;
% Temporary dir
tmpdir = conf.tmpdir;


%% ===== Simulation =====================================================
t = round(linspace(-20,900,10*25));
% Generate a random number string for the tmp files
rn = sprintf('%04.0f',10000*rand);
% Disable the empty sound field warning
warning('off','SFS:check_sound_field');
conf.plot.useplot = false;
conf.plot.usenormalisation = false;
% Simulate the time by different phase values
for ii = 1:length(t)-1
    % Calculate sound field for the given phase
    [p,x,y,z,x0] = sound_field_imp_wfs(X,Y,Z,xs,src,t(ii),conf);

    % === Save temporary data ===
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end
    conf.plot.file = sprintf('%s/%s_%04.0f.png',tmpdir,rn,ii);
    %conf.plot.usedb = 1;
    plot_sound_field(2*p,x,y,z,x0,conf);
end
% Enable the empty sound field warning
warning('on','SFS:check_sound_field');

%% ===== Create movie ====================================================
conf.plot.useplot = useplot;
generate_movie(outfile,tmpdir,rn);

% Clean up tmp files
delete([tmpdir,'/',rn,'*.png']);


%% ===== Show movie ======================================================
if useplot
    status = system('which mplayer');
    if status
        error('%s: mplayer is needed to show this movie.',upper(mfilename));
    else
        cmd = sprintf('mplayer %s -loop 0',outfile);
        system(cmd);
    end
end
