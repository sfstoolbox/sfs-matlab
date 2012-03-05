function [s] =  bandpass(s,conf)
%BANDPASS filters a signal by a bandpass
%
%   Usage: s = bandpass(s,conf)
%          
%
%
%   Input parameters:
%       s       - input signal (vector)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       s    - delayed signal

% AUTHOR: Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$



%% ===== Checking of input  parameters ==================================


%% ===== Configuration ==================================================
N = 256;
%N = 4*4096;

fs = conf.fs;

%% ===== Computation =====================================================    

% design bandpass filter
Hf = [0 2*10/fs 2*20/fs 2*18000/fs 2*20000/fs 1]; 
Hm = [0 0 1 1 0 0];
b = fir2(N,Hf,Hm);

% filter signal
s = conv(s,b);

% compensate for delay & truncate result
s = s(N/2:end-(N/2)-1);