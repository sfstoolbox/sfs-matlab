function generate_movie(outfile,directory,pattern)
%GENERATE_MOVIE generates a movie from given png files
%
%   Usage: generate_movie(outfile,[[directory],pattern])
%
%   Input parameters:
%       outfile     - file name for the movie
%       directory   - dir containing the png files
%       pattern     - starting pattern of the png files to use
%
%   GENERATE_MOVIE(outfile,directory,pattern) generates a movie and stores it
%   in outfile. To produce the movie all png files starting with pattern in the
%   given directory are processed. If no pattern is given, all png files are
%   used.
%
%   See also: plot_sound_field, movie_sound_field_mono_wfs

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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
nargmin = 1;
nargmax = 3;
narginchk(nargmin,nargmax);
if ~exist('directory','var')
    directory = pwd;
end
if ~exist('pattern','var')
    pattern = '';
end
isargchar(outfile,pattern);
isargdir(directory);


%% ===== Movie ===========================================================
% Generate a movie from the png files with MEncoder
status1 = system('which mencoder > /dev/null');
status2 = system('which avconv > /dev/null');
if status1 & status2
    error('%s: mencoder or avconv is needed to generate the movie.', ...
        upper(mfilename));
else
    if ~status1
        cmd = sprintf(['mencoder mf://%s/%s*.png -mf fps=25:type=png', ...
                       ' -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell ', ...
                       '-oac copy -o %s'],...
                      directory,pattern,outfile);
    else
        cmd = sprintf(['avconv -f image2 -r 25 -i %s/%s_%%04d.png ', ...
                       '-c:v libx264 -crf 19 %s'], ...
                      directory,pattern,outfile);
    end
    system(cmd);
end
