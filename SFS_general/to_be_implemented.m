function to_be_implemented()
%TO_BE_IMPLEMENTED indicates that the code had to be implemented
%   Usage: to_be_implemented(mfile)
%          to_be_implemented()
%
%   Input parameters:
%       mfile   - string containing the name of the calling m file.
%                 NOTE: this variable is already avaiable in all m files as
%                 mfilename!
%
%   TO_BE_IMPLEMENTED(mfile) results in an error that indicates the desired code
%   functionality has to be implemented yet.

% AUTHOR: Hagen Wierstorf

%% ===== Checking of input parameters ====================================
%error(nargchk(0,1,nargin)
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
    error('%s: This functionality has to be implemented yet!',upper(mfile));
else
    error('This functionality has to be implemented yet!');
end
