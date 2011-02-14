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

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
% Disable the ordering warning, becuase we wanted to reorder the entries
warning('off','SFS:irs_fields_order');
check_irs(irs);
warning('on','SFS:irs_fields_order');


%% ===== Computation =====================================================
% Get the reference implementation of the irs format and reorder the fields of
% the given irs according to the reference implementation.
ref_irs = new_irs();
% Get fields
ref_fields = fieldnames(ref_irs);
% Get the fields for the given irs
fields = fieldnames(irs);
% If the number of fields is identical sort directly
if length(ref_fields)==length(fields)
    irs = orderfields(irs,ref_fields);
else
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
