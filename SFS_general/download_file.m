function success = download_file(url,outfile)
%DOWNLOAD_FILE downloads a file to your computer
%
%   Usage: success = download_file(url,outfile)
%
%   Input parameters:
%       url     - url to download
%       outfile - path to store the file
%
%   Output parameters:
%       success  - 0 or 1
%
%   DOWNLOAD_FILE(url,file) downloads the given url and stores it at outfile.
%   If outfile contains directories that do not exist yet, they will be created.
%
%   See also: get_spherical_grid

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


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargchar(url,outfile)


%% ===== Main ============================================================
% Download if file is not present
if ~exist(outfile,'file')
    % Replace '\' with '/'
    outfile = strrep(outfile,'\','/');
    % Create dir
    % NOTE: newer versions of Matlab can do the following with the strsplit
    % function
    %dirs = strsplit(outfile,'/');
    dirs = regexp(outfile,'/','split');
    dir_path = '';
    for ii=1:length(dirs)-1
        if ii==1 && iswindows
            dir_path = [dir_path dirs{ii}];
        else
            dir_path = [dir_path '/' dirs{ii}];
        end
        [~,~] = mkdir(dir_path);
    end
    warning('Downloading file %s',url);
    [~,success] = urlwrite(url,outfile);
else
    error('%s: file exist.',upper(mfilename));
end
