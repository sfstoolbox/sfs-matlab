function inoutsig = gaindb(inoutsig,gn)
%GAINDB  Increase/decrease level of signal
%   Usage:  outsig = gaindb(insig,gn);
%
%   Input parameters:
%       inoutsig    - signal for which the level should be changed
%       gn          - gain of the level (dB)
%
%   Output parameters:
%       inoutsig    - given signal with new level
%
%   GAINDB(insig,gn) changes the level of the signal by gn dB.
%
%   See also: rms, db, setleveldb
%

%   AUTHOR: Hagen Wierstorf, Peter L. Soendergaard


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargnumeric(inoutsig)
isargscalar(gn)


%% ===== Computation =====================================================
inoutsig = inoutsig*10^(gn/20);
