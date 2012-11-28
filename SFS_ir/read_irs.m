function irs = read_irs(irsfile)
%READ_IRS Read a HRIR/BRIR dataset
%
%   Usage: irs = read_irs(irsfile)
%
%   Input parameters:
%       irsfile - filename of irs mat file
%
%   Output paramteres:
%       irs   - irs struct. For details on the containing fields have a look at
%               the IR_format.txt file.
%
%   READ_IRS(irsfile) loads a IR dataset as a struct containing the format
%   specific fields. For a description of the mat format for the IR datasets,
%   see IR_format.txt.
%
%   see also: get_ir, intpol_ir, dummy_irs, new_irs, brs_wfs_25d,
%       brs_point_source

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargfile(irsfile);


%% ===== Read IR files ================================================
% Load the mat file
load(irsfile);
% Check irs format
check_irs(irs);
