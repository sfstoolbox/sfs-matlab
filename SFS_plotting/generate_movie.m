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
