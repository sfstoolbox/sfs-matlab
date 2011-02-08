function sigdb = db(sig)
%DB applies 20*log10() to the given signal
%   Usage: sigdb = db(sig);
%
%   Input options:
%       sig     - signal to calculate the level
%
%   Output options:
%       level   - level of given signal (dB). It has the same size as the input
%                 signal
%
%   DB(sig) returns the level in dB for the given input values by applying a
%   simple 20*log10(sig).
%
%   see also: rmsdb

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(sig);


%% ===== Computation =====================================================
sigdb = 20*log10(sig);
