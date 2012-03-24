function label = print_label(dim,unit,conf)
%PRINT_LABEL retruns a suitable axis label
%
%   Usage: label = print_label(dim,[unit,[conf]])
%
%   Input parameters:
%       dim     - name of label dimension
%       unit    - name of label unit, default: no unit
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Ouput parameters:
%       label   - label to put on an axis of a plot
%
%   PRINT_LABEL(DIM,UNIT) generates a label with the given axis dimension and
%   unit. The formatting is depending on your plotting style, adding $$ for
%   LaTeX labels.
%
%   see also:

% AUTHOR: Hagen Wierstorf
% $LastChangedDate: $
% $LastChangedRevision: $
% $LastChangedBy: $


%% ===== Checking of input parameter =====================================
nargmin = 1;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
if nargin==nargmax-1
    if isstruct(unit)
        conf = unit;
        unit = '';
    else
        conf = SFS_config;
    end
elseif nargin==nargmax-2
    unit = '';
    conf = SFS_config;
end
isargchar(dim,unit);
isargstruct(conf);


%% ===== Configuration ===================================================
p.mode = conf.plot.mode;


%% ===== Main ============================================================
if strcmp(p.mode,'paper') | strcmp(p.mode,'talk')
    if length(unit)>0
        label = sprintf('$%s$~(%s)',dim,unit);
    else
        label = sprintf('$%s$',dim);
    end
else
    if length(unit)>0
        label = sprintf('%s (%s)',dim,unit);
    else
        label = sprintf('%s',dim);
    end
end
