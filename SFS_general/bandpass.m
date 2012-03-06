function [s] =  bandpass(s,conf)
%BANDPASS filters a signal by a bandpass
%
%   Usage: s = bandpass(s,conf)
%          
%   Input parameters:
%       s    - input signal (vector)
%       conf - optional struct containing configuration variables (see
%              SFS_config for default values)
%
%   Output parameters:
%       s    - delayed signal
%
%   BANDPASS(sig) filters the given signal with a bandpass filter
%   with ...

% AUTHOR: Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargvector(s);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
fs = conf.fs;
N = 256;
%N = 4*4096;


%% ===== Computation =====================================================    
% design bandpass filter
Hf = [0 2*10/fs 2*20/fs 2*18000/fs 2*20000/fs 1]; 
Hm = [0 0 1 1 0 0];
b = fir2(N,Hf,Hm);
% filter signal
s = conv(s,b);
% compensate for delay & truncate result
s = s(N/2:end-(N/2)-1);
