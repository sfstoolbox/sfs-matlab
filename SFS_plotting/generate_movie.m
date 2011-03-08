function generate_movie(outfile,directory,pattern)
%GENERATE_MOVIE generates a movie from given png files
%   Usage: generate_movie(outfile,directory,pattern)
%          generate_movie(outfile,directory)
%          generate_movie(outfile)
%
%   Input parameters:
%       outfile     - file name for the movie
%       directory   - dir containing the png files
%       pattern     - starting pattern of the png files to use
%
%   GENERATE_MOVIE(outfile,directory,pattern) generates a movie and stores it
%   in outfile. To produce the movie all png files starting with pattern in the
%   given directory are processed. If no pattern is given, all png files are
%   used. framerate?
%
%   see also: plot_wavefield, movie_wave_field_mono_wfs_25d 
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
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
[status,encoder] = system('which mencoder');
if status
    error('%s: mencoder is needed to generate the movie.',upper(mfilename));
else
    cmd = sprintf(['mencoder mf://%s/%s*.png -mf fps=25:type=png', ...
        ' -ovc copy -oac copy -o %s'],directory,pattern,outfile);
    system(cmd);
end
