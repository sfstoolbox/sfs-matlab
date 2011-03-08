function P = norm_wave_field(P,x,y,conf)
%NORM_WAVE_FIELD normalizes the wave field to 1 at xref,yref
%   Usage: P = norm_wave_field(P,x,y,yref)
%
%   Input options:
%       P       - wave field
%       x,y     - vectors conatining the x- and y-axis values
%       yref    - yref coordinate for normalization
%
%   Output options:
%       P       - normalized wave field
%
%   NORM_WAVE_FIELD(P,x,y,yref) normalizes the given wave field P to 1 at
%   the position x/2,yref.
%
%   see also: wave_field_mono_wfs_25d

% AUTHOR: Hagen Wierstorf


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
xref = conf.xref;
yref = conf.yref;


%% ===== Computation =====================================================
% Use the half of the x axis and yref
[a,xidx] = find(x>xref,1);
[a,yidx] = find(y>yref,1);
if isempty(xidx)
    error('%s: your used xref is out of your X boundaries',upper(mfilename));
end
if isempty(yidx)
    error('%s: your used yref is out of your Y boundaries',upper(mfilename));
end
% Scale signal to 1
P = 1*P/abs(P(yidx,xidx));
