function movie_sound_field_mono_wfs(X,Y,Z,xs,src,f,outfile,conf)
%MOVIE_SOUND_FIELD_MONO_WFS_25D generates movie a WFS sound field
%
%   Usage: movie_sound_field_mono_wfs_25d(X,Y,Z,xs,src,f,outfile,[conf])
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
%       f           - monochromatic frequency / Hz
%       outfile     - name for the movie file
%       conf        - optional configuration struct (see SFS_config)
%
%   MOVIE_SOUND_FIELD_MONO_WFS(X,Y,Z,xs,src,f,L,outfile,conf) generates a
%   movie of simulations of a sound field of the given source positioned at xs
%   using a WFS driving function in the temporal domain with different phase.
%
%   see also: sound_field_mono_wfs, plot_sound_field

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
nargmin = 7;
nargmax = 8;
narginchk(nargmin,nargmax);
isargvector(X,Y,Z);
isargxs(xs);
isargpositivescalar(f);
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
phase = linspace(2*pi,0,25);
% Generate a random number string for the tmp files
rn = sprintf('%04.0f',10000*rand);
% Simulate the time by different phase values
for ii = 1:length(phase)-1
    conf.phase = phase(ii);
    conf.plot.useplot = 0;
    % Calculate sound field for the given phase
    [P,x,y,z,x0] = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);

    % === Save temporary data ===
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end
    conf.plot.file = sprintf('%s/%s_%i.png',tmpdir,rn,ii+10);
    plot_sound_field(P,x,y,z,x0,conf);
end


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
