function movie_sound_field_imp_wfs(X,Y,Z,xs,src,outfile,conf)
%MOVIE_SOUND_FIELD_IMP_WFS generates movie a WFS sound field
%
%   Usage: movie_sound_field_imp_wfs_25d(X,Y,Z,xs,src,outfile,[conf])
%
%   Input parameters:
%       X           - x-axis / m; single value or [xmin,xmax] 
%       Y           - y-axis / m; single value or [ymin,ymax]
%       Z           - z-axis / m; single value or [zmin,zmax]
%       xs          - position of point source / m
%       src         - sourcetype of the virtual source:
%                         'pw' - plane wave (xs, ys are the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       outfile     - name for the movie file
%       conf        - optional configuration struct (see SFS_config)
%
%   MOVIE_SOUND_FIELD_IMP_WFS(X,Y,Z,xs,src,outfile,conf) generates a movie of
%   simulations of a sound field of the given source positioned at xs
%   using a WFS driving function in the temporal domain with different phase.
%
%   see also: sound_field_imp_wfs, plot_sound_field, generate_movie

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


%% ===== Checking of input  parameters ==================================
nargmin = 6;
nargmax = 7;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z);
isargxs(xs);
isargchar(src,outfile);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


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
conf.plot.useplot = 0;
conf.usenormalisation = 0;
% Simulate the time by different phase values
for ii = 1:length(t)-1
    % Calculate sound field for the given phase
    [p,x,y,z,x0,win] = sound_field_imp_wfs(X,Y,Z,xs,src,t(ii),conf);

    % === Save temporary data ===
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end
    conf.plot.file = sprintf('%s/%s_%04.0f.png',tmpdir,rn,ii);
    %conf.plot.usedb = 1;
    plot_sound_field(2*p,x,y,z,x0,win,conf);
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
