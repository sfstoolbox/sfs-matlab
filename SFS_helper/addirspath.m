function addirspath(conf)
%ADDIRSPATH adds directories containing irs files to the path
%
%   Usage: addirspath(conf)
%
%   Input parameters:
%       conf    - optional configuration struct (see SFS_config)
%
%   ADDIRSPATH(conf) adds the directories specified in conf.ir.dir to your
%   search path for easy loading with read_irs(). If you have more than one
%   directory they have to be seperated by ':'.

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


% ===== Checking of input parameters ====================================
nargmin = 0;
nargmax = 1;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    conf = SFS_config;
end
isargstruct(conf);


%% ===== Configuration ==================================================
% NOTE: newer versions of Matlab can do the following with the strsplit
% function
%dirs = strsplit(conf.ir.path,':');
dirs = regexp(conf.ir.path,':','split')


%% ===== Adding pathes ==================================================
for ii = 1:length(dirs)
    if exist(dirs{ii},'dir')
        path = genpath(dirs{ii});
        addpath(path,'-end');
    else
        warning('%s: %s does not exist.',upper(mfilename),dirs{ii});
    end
end
