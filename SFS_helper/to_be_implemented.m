function to_be_implemented(mfile)
%TO_BE_IMPLEMENTED indicates that the code had to be implemented
%   Usage: to_be_implemented(mfilename)
%          to_be_implemented()
%
%   Input parameters:
%       mfilename   - string containing the name of the calling m file.
%                   NOTE: this variable is already avaiable in all m files as
%                   mfilename!
%
%   TO_BE_IMPLEMENTED(mfilename) results in an error that indicates the desired 
%   code functionality has to be implemented yet.

% AUTHOR: Hagen Wierstorf

%% ===== Checking of input parameters ====================================
nargmin = 0;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin);


%% ===== Main ============================================================
if exist('mfile','var')
    isargchar(mfile)
    error('%s: This functionality has to be implemented yet!',upper(mfile));
else
    error('This functionality has to be implemented yet!');
end
