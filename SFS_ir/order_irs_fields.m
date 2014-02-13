function irs = order_irs_fields(irs)
%ORDER_IRS_FIELDS orders the fields in irs according to the new_irs reference
%
%   Usage: irs = order_irs_fields(irs)
%
%   Input options:
%       irs -   irs struct
%
%   Output options
%       irs - irs struct with the fields in the order as in new_irs()
%
%   ORDER_IRS_FIELDS(irs) reorders the fields in the given irs according to the
%   order given in the new_irs function, which is our reference implementation
%   for the irs format. If the given irs struct contains more fields, the
%   additional fields will be appended at the end of the struct in the order
%   given in the original struct.
%
%   See also: new_irs, check_irs, IR_format.txt

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
% Disable the ordering warning, becuase we wanted to reorder the entries
warning('off','SFS:irs_fields_order');
check_irs(irs);
warning('on','SFS:irs_fields_order');


%% ===== Computation =====================================================
% Get the reference implementation of the irs format and reorder the fields of
% the given irs according to the reference implementation.
[ref_irs,opt_fields] = new_irs();
% Get fields
ref_fields = fieldnames(ref_irs);
% Get the fields for the given irs
fields = fieldnames(irs);
% If the number of fields is identical sort directly
if length(ref_fields)==length(fields)
    irs = orderfields(irs,ref_fields);
else
    % Remove unneeded optional fields from the reference
    idx = [];
    for ii = 1:length(ref_fields)
        if ~isfield(irs,ref_fields{ii}) && strcmp(opt_fields,ref_fields{ii})
            idx = [idx ii];
        end
    end
    ref_fields(idx) = [];
    % Get the indices of the positions of the ref_fields in the given irs
    idx = [];
    for ii = 1:length(ref_fields)
        for jj = 1:length(fields)
            if strcmp(fields{jj},ref_fields{ii})
                idx = [idx,jj];
            end
        end
    end
    % Get the indices for additional fields
    for jj = 1:length(fields)
        if ~any(idx==jj)
            idx = [idx,jj];
        end
    end
    sorted_fields = fields(idx);
    irs = orderfields(irs,sorted_fields);
end
