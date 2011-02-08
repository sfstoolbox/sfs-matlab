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
%   ORDER_IRS_FIELDS() reorders the fields in the given irs according to the
%   order given in the new_irs function, which is our reference implementation
%   for the irs format.
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
irs = orderfields(irs,ref_irs);
