function sig = norm_signal(sig)
%NORM_SIGNAL normalizes the signal to 1-eps
%
%   Usage: sig = norm_signal(sig)
%
%   Input parameters:
%       sig     - audio signal
%
%   Output parameters:
%       sig     - audio signal normalized to -1<sig<1
%
%   NORM_SIGNAL(sig) normalizes the amplitude of the given signal to the range
%   of ]-1:1[.
%
%   see also: auralize_ir
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(sig);


%% ===== Computation =====================================================
% Scaling the signal to -1<sig<1
sig = sig / (max(abs(sig(:)))+eps);
