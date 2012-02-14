function P = norm_wave_field(P,x,y,conf)
%NORM_WAVE_FIELD normalizes the wave field to 1 at xref,yref
%   Usage: P = norm_wave_field(P,x,y,conf)
%
%   Input options:
%       P       - wave field
%       x,y     - vectors conatining the x- and y-axis values
%       conf    - optional configuration struct
%
%   Output options:
%       P       - normalized wave field
%
%   NORM_WAVE_FIELD(P,x,y,yref) normalizes the given wave field P to 1 at
%   the position conf.xref,conf.yref.
%
%   see also: wave_field_mono_wfs_25d

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 3;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(P);
isargvector(x,y);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
xref = position_vector(conf.xref);


%% ===== Computation =====================================================
% Use the half of the x axis and xref
[a,xidx] = find(x>xref(1),1);
[a,yidx] = find(y>xref(2),1);
if isempty(xidx) || abs(x(xidx)-xref(1))>0.1
    error('%s: your used conf.xref is out of your X boundaries',upper(mfilename));
end
if isempty(yidx) || abs(y(yidx)-xref(2))>0.1
    error('%s: your used conf.xref is out of your Y boundaries',upper(mfilename));
end
% Scale signal to 1
P = 1*P/abs(P(yidx,xidx));
