function ir = shorten_ir(ir,nsamples)
%SHORTEN_IR shortens a IR
%   Usage: ir = shorten_ir(ir,nsamples)
%
%   Input parameters:
%       ir          - IR signal with length x channels
%       nsamples    - length of the target IR
%
%   Output paramteres:
%       ir          - IR signal with nsamples x n
%
%   SHORTEN_HRIR(ir,nsamples) shortens a given IR to the given number of samples
%   nsamples and applying a 5% long hanning window.
%
%   see also: SFS_config, read_irs, intpol_ir
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargpositivescalar(nsamples);


%% ===== Computation ====================================================

% Window IR
win = hanningwin(0,ceil(0.05*nsamples),nsamples);

ir = ir(1:nsamples,:) .* repmat(win,[1 size(ir,2)]);
