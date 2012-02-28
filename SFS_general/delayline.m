function [s] =  delayline(s,dt,weight,hf)
%DELAYLINE implements a (fractional) delay line
%
%   Usage: s = delayline(s,dt,weight,hf)
%          
%
%
%   Input parameters:
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       s    - delayed signal

% AUTHOR: Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


% apply integer delay first
idt=floor(dt);
s = [zeros(1,idt) weight*s(1:end-idt)];


if(nargin>3)
    if(~isempty(hf))
    hf.FracDelay = dt-idt;
    s = filter(hf,s);
    end
end