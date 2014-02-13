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
%   see also: get_spherical_grid

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


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargchar(url,outfile)


%% ===== Main ============================================================
% download if file is not present
if ~exist(outfile,'file')
    % replace '\' with '/'
    outfile = strrep(outfile,'\','/');
    % create dir
    % NOTE: newer versions of Matlab can do the following with the strsplit
    % function
    %dirs = strsplit(outfile,'/');
    dirs = regexp(outfile,'/','split');
    dir_path = [];
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
