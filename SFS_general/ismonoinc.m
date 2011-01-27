function bool = ismonoinc(x)
%ISMONOINC Checks if the data of x are monotonic increasing
%   Usage: bool = ismonoinc(x)
%
%   Input parameter:
%       x       - data vector
%
%   Output parameter:
%       bool    - boolean value
%
%   ISMONOTONIC(x) checks if the data in vector x are monotonic increasing.
%
%   see also: hrir_intpol
%

% AUTHOR: Hagen Wierstorf

%% ===== Checking of input  parameters ==================================

if nargchk(1,1,nargin)
    error('Wrong number of args. Usage: bool = ismonoinc(x)');
end

if ~isnumeric(x) || ~isvector(x)
    error('%s: x has to be a vector!',upper(mfilename));
end


%% ===== Computation ====================================================
bool = true;
for i = 2:length(x)
    if x(i)<x(i-1)
        bool = false;
        return;
    end
end


