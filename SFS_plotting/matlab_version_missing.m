function matlab_version_missing(mfile)
%MATLAB_VERSION_MISSING indicates a function that is not implemented yet
%   Usage: matlab_version_missing(mfilename)
%          matlab_version_missing()
%
%   Input parameters
%       mfile   - string containing the name of the calling m file.
%                 NOTE: this variable is already avaiable in all m files as
%                 mfilename!
%
%   MATLAB_VERSION_MISSING(mfile) results in an error that indicates that the
%   desired plot functionality has to be implemented in Matlab. This is the 
%   case for a some heavier plotting functions that are first implemented using 
%   Gnuplot.

% AUTHOR: Hagen Wierstorf

%% ===== Checking of input parameters ====================================
error(nargchk(0,1,nargin)
if exist('mfile','var') 
    if ~ischar(mfile)
        error('%s: mfile has to be a string.',upper(mfilename));
    end
    filename = true;
else
    filename = false
end


%% ===== Main ============================================================
if(filename)
    error(['%s: The Matlab version of the movie function has to be ', ...
        'implemented yet. You can use the Gnuplot version by setting ', ...
        'conf.gnuplot = 1. If you need the Matlab version, please feel ', ...
        'free to implement it ;)'],upper(mfile));
else
    error(['The Matlab version of the movie function has to be ', ...
        'implemented yet. You can use the Gnuplot version by setting ', ...
        'conf.gnuplot = 1. If you need the Matlab version, please feel ', ...
        'free to implement it ;)']);
end
