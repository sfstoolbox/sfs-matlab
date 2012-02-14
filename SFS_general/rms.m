function y = rms(insig,options)
%RMS returns the root mean square of the signal
%
%   Usage: y = rms(insig);
%          y = rms(insig,'ac');
%
%   Input parameters:
%       insig       - signal for which to calculate the RMS value
%       options     - if 'ac' is given, only the AC component of the signal
%                     is used for RMS calculation
%
%   Output parameters:
%       y           - RMS value of insig
%
%   RMS(insig) computes the RMS (Root Mean Square) value of a finite 
%   sampled signal sampled at a uniform sampling rate.
%
%   RMS(x,'ac') does the same, but considers only the AC component of the
%   signal (i.e. the mean is removed).
%
%   The RMS value of a signal x of length N is computed by
%
%                        1  N
%     rms(insig) = sqrt( - sum insig(n)^2 )
%                        N n=1
%

% AUTHOR : Hagen Wierstorf, Peter L. Soendergaard
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargvector(insig);
if (nargin==1) || (~ischar(options))
  options='';
end


%% ===== Computation =====================================================
% It is better to use 'norm' instead of explicitly summing the squares, as
% norm (hopefully) attempts to avoid numerical overflow.
switch(lower(options))
    case 'ac'
        y = norm(insig-mean(insig))/sqrt(length(insig));
    otherwise
        y = norm(insig)/sqrt(length(insig));
end
