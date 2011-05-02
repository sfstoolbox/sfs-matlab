function h0 = firintpolc(n,omi,Hw,g)
%firintpolc: Design of FIR filter for a ideal complex frequency response
%
%   Syntax:
%   h0 = firintpolc(n,omi,Hw,g)
%
%   Description:
% h0 = firintpolc(n,omi,Hw,g) calculates the impulse response of a FIR 
%   filter of order n with real coefficients. The frequency response of 
%   that filter approximates the ideal function Hw(exp (j*omi*pi)). 
%   The ideal function is defined at the normalized frequency points 
%   omi(k) = Omega(k)/pi; omi=[0 1). G(k) is an optional weighting function
%   for the design. By default G(k)is set to G(k)=1.
%   The vectors omi, Hw and g have to be of the same length. 
%
%   Algorithm: Design by interpolation
%
%   Remark:
%
%   References: [DSV2, chap. 2.11.3]
%
%   See also  DSV-Lib: firminGrpL2c, firRGcTschebyLp
%             Matlab:  mldivide, matrix left division

%% ---------------------------------------------------------------
%   DSV2-Bibliothek zu H.W. Schuessler: Digitale Signalverarbeitung 2,
%                      Springer (2010)
%   Authors: Hans W. Schuessler and Guenter F. Dehner
%   Version: 1.0        Release DSV 2010/1
%   Copyright 2009 by authors - not released for commercial use
%%  ---------------------------------------------------------------
error(nargchk(3,4,nargin, 'struct'))
%
m = length(omi); k = 0:n; omi = omi(:); Hw = Hw(:);
if nargin == 4, g = g(:); else g = ones(m,1); end
if omi(1) == 0,
   omig = [-flipud(omi); omi(2:m)];
   Hwg = [flipud(conj(Hw)); Hw(2:m)];
   gg = [flipud(g); g(2:m)];
else
   omig = [-flipud(omi);omi];
   Hwg = [flipud(conj(Hw));Hw];
   gg = [flipud(g); g];
end
A = exp(-1i*omig*pi*k);
G = diag(gg);
h0 = real(G*A\(G*Hwg)); 
h0=h0(:).';
%% [EOF] - firintpolc

