function bool = ismonoinc(x)
%ISMONOINC checks if the data of x are monotonic increasing
%
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
%   see also: intpol_ir
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
isargvector(x);


%% ===== Computation ====================================================
bool = true;
for i = 2:length(x)
    if x(i)<x(i-1)
        bool = false;
        return;
    end
end
