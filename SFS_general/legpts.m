function [x w v] = legpts(n,int,meth)
%LEGPTS  Legendre points and Gauss Quadrature Weights.
%  LEGPTS(N) returns N Legendre points X in (-1,1).
%
%  [X,W] = LEGPTS(N) returns also a row vector W of weights for Gauss quadrature.
%
%  LEGPTS(N,D) scales the nodes and weights for the domain D. D can be
%  either a domain object or a vector with two components. If the interval
%  is infinite, the map is chosen to be the default 'unbounded map' with
%  mappref('parinf') = [1 0] and mappref('adaptinf') = 0.
%
%  [X,W,V] = LEGPTS(N) returns additionally a column vector V of weights in
%  the barycentric formula corresponding to the points X.
%
%  [X,W] = LEGPTS(N,METHOD) allows the user to select which method to use.
%    METHOD = 'FASTSMALL' uses the recurrence relation for the Legendre 
%       polynomials and their derivatives to perform Newton iteration 
%       on the WKB approximation to the roots.
%    METHOD = 'FAST' uses the Glaser-Liu-Rokhlin fast algorithm, which
%       is much faster for large N. 
%    METHOD = 'GW' will use the traditional Golub-Welsch eigenvalue method, 
%       which is maintained mostly for historical reasons.
%       By default LEGPTS uses 'FASTSMALL' when N<=256 and FAST when N>256
%
%  See also chebpts and jacpts.

%  Copyright 2011 by The University of Oxford and The Chebfun Developers. 
%  See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%      * Neither the name of the University of Oxford nor the names of its 
%        contributors may be used to endorse or promote products derived from 
%        this software without specific prior written permission.
%  
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR 
%  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%  'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
%  'FAST' by Nick Hale, April 2009 - algorithm adapted from [2].
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature
%       rules", Math. Comp. 23:221-230, 1969, 
%   [2] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the 
%       calculation of the roots of special functions", SIAM Journal  
%       on Scientific Computing", 29(4):1420-1438:, 2007.

% Defaults
interval = [-1,1];
method = 'default';

if n < 0
    error('CHEBFUN:legpts:n','First input should be a positive number.');
end

% Return empty vector if n == 0
if n == 0
    x = []; w = []; v = []; return
end

% Check the inputs
if nargin > 1
    if nargin == 3
        interval = int; method = meth;
    elseif nargin == 2
        if ischar(int), method = int; else interval = int; end
    end
    if ~any(strcmpi(method,{'default','GW','fast','fastsmall'}))
        error('CHEBFUN:legpts:inputs',['Unrecognised input string.', method]); 
    end
    if isa(interval,'domain')
        interval = interval.endsandbreaks;
    end
    if numel(interval) > 2,
        warning('CHEBFUN:legpts:domain',...
            'Piecewise intervals not supported and will be ignored.');
        interval = interval([1 end]);
    end
end

% Decide to use GW or FAST
if n == 1
% Trivial case when N = 1    
    x = 0; w = 2; v = 1;
elseif strcmpi(method,'GW')
% GW, see [1]
   beta = .5./sqrt(1-(2*(1:n-1)).^(-2)); % 3-term recurrence coeffs
   T = diag(beta,1) + diag(beta,-1);     % Jacobi matrix
   [V,D] = eig(T);                       % Eigenvalue decomposition
   x = diag(D); [x,i] = sort(x);         % Legendre points
   w = 2*V(1,i).^2;                      % Quadrature weights
   v = sqrt(1-x.^2).*abs(V(1,i))';       % Barycentric weights
   v = v./max(v);
   
   % Enforce symmetry
   ii = 1:floor(n/2);  x = x(ii);  w = w(ii); 
   vmid = v(floor(n/2)+1); v = v(ii);
   if mod(n,2)
        x = [x ; 0 ; -x(end:-1:1)];  w = [w  2-sum(2*w) w(end:-1:1)];
        v = [v ; vmid ; v(end:-1:1)];
   else
        x = [x ; -x(end:-1:1)];      w = [w w(end:-1:1)];      
        v = [v ; v(end:-1:1)];
   end
   v(2:2:n) = -v(2:2:end);
elseif (n < 256 && ~strcmpi(method,'fast')) || strcmpi(method,'fastsmall')
% Fastsmall    
   [x ders] = fastsmall(n);              % Nodes and P_n'(x)
   w = 2./((1-x.^2).*ders.^2)';          % Quadrature weights
   v = 1./ders; v = v./max(abs(v));      % Barycentric weights  
   if ~mod(n,2), ii = (floor(n/2)+1):n; v(ii) = -v(ii);   end
else
% Fast, see [2]
   [x ders] = alg0_Leg(n);               % Nodes and P_n'(x)
   w = 2./((1-x.^2).*ders.^2)';          % Quadrature weights
   v = 1./ders; v = v./max(abs(v));      % Barycentric weights
   if ~mod(n,2), ii = (floor(n/2)+1):n;  v(ii) = -v(ii);   end
end

% Normalise so that sum(w) = 2
w = (2/sum(w))*w;           

% Rescale to arbitrary interval
if ~all(interval == [-1 1])
    if ~any(isinf(interval))
    % Finite interval
        dab = diff(interval);
        x = (x+1)/2*dab + interval(1);
        w = dab*w/2;
    else
    % Infinite interval
        m = maps(fun,{'unbounded'},interval); % use default map
        if nargout > 1
            w = w.*m.der(x.');
        end
        x = m.for(x);
        x([1 end]) = interval([1 end]);
    end
end

% -------------------- Routines for FAST_SMALL algorithm ------------------

function [x PP] = fastsmall(n)

% Asymptotic formula (WKB) - only positive x.
if mod(n,2), s = 1; else s = 0; end 
k = (n+s)/2:-1:1; theta = pi*(4*k-1)/(4*n+2);
x = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);

% Initialise
Pm2 = 1; Pm1 = x;  PPm2 = 0; PPm1 = 1;
dx = inf; l = 0;

% Loop until convergence
while norm(dx,inf) > eps && l < 10
    l = l + 1;
    for k = 1:n-1, 
        P = ((2*k+1)*Pm1.*x-k*Pm2)/(k+1);           Pm2 = Pm1; Pm1 = P; 
        PP = ((2*k+1)*(Pm2+x.*PPm1)-k*PPm2)/(k+1);  PPm2 = PPm1; PPm1 = PP;  
    end
    dx = -P./PP; x = x + dx;    
    Pm2 = 1; Pm1 = x; PPm2 = 0; PPm1 = 1;
end

% Once more for derivatives
for k = 1:n-1, 
    P = ((2*k+1)*Pm1.*x-k*Pm2)/(k+1);           Pm2 = Pm1; Pm1 = P; 
    PP = ((2*k+1)*(Pm2+x.*PPm1)-k*PPm2)/(k+1);  PPm2 = PPm1; PPm1 = PP;  
end

% Reflect for negative values
x = [-x(end:-1:1+s) x].';
PP = [PP(end:-1:1+s) PP].';


% -------------------- Routines for FAST algorithm ------------------------

function [roots ders] = alg0_Leg(n) % Driver for 'Fast'.

% Compute coefficients of P_m(0), m = 0,..,N via recurrence relation.
Pm2 = 0; Pm1 = 1; 
for k = 0:n-1, P = -k*Pm2/(k+1); Pm2 = Pm1; Pm1 = P; end

% Get the first roots and derivative values to initialise.
roots = zeros(n,1); ders = zeros(n,1);                      % Allocate storage
if mod(n,2)                                                 % n is odd
    roots((n-1)/2) = 0;                                     % Zero is a root
    ders((n+1)/2) = n*Pm2;                                  % P'(0)    
else                                                        % n is even
    [roots(n/2+1) ders(n/2+1)] = alg2_Leg(P,n);             % Find first root
end       

[roots ders] = alg1_Leg(roots,ders);          % Other roots and derivatives

% -------------------------------------------------------------------------

function [roots ders] = alg1_Leg(roots,ders)  % Main algorithm for 'Fast'
n = length(roots);
if mod(n,2), N = (n-1)/2; s = 1; else N = n/2; s = 0; end   

% Approximate roots via asymptotic formula.
k = (n-2+s)/2:-1:1; theta = pi*(4*k-1)/(4*n+2);
roots(((n+4-s)/2):end) = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);
x = roots(N+1);

% Number of terms in Taylor expansion.
m = 30;

% Storage
hh1 = ones(m+1,1); zz = zeros(m,1); u = zeros(1,m+1); up = zeros(1,m+1);

% Loop over all the roots we want to find (using symmetry).
for j = N+1:n-1
    % Distance to initial approx for next root (from asymptotic foruma).
    h = roots(j+1) - x;

    % Recurrence Taylor coefficients (scaled & incl factorial terms).
    M = 1/h;                           % Scaling
    c1 = 2*x/M; c2 = 1./(1-x^2);       % Some constants
    % Note, terms are flipped for more accuracy in inner product calculation.
    u([m+1 m]) = [0 ders(j)/M];  up(m+1) = u(m);
    for k = 0:m-2
        up(m-k) = (c1*(k+1)*u(m-k)+(k-n*(n+1)/(k+1))*u(m-k+1)/M^2)*c2;
        u(m-(k+1)) = up(m-k)/(k+2);
    end
    up(1) = 0;  

    % Newton iteration
    hh = hh1; step = inf;  l = 0; 
    while (abs(step) > eps) && (l < 10)
        l = l + 1;
        step = (u*hh)/(up*hh)/M;
        h = h - step;        
        Mhzz = (M*h)+zz;
        hh = [1;cumprod(Mhzz)];     % Powers of h (This is the fastest way!)
        hh = hh(end:-1:1);          % Flip for more accuracy in inner product 
    end
      
    % Update
    x = x + h;
    roots(j+1) = x;
    ders(j+1) = M*(up*hh);  
       
end

% Nodes are symmetric.
roots(1:N+s) = -roots(n:-1:N+1);
ders(1:N+s) = ders(n:-1:N+1);

% -------------------------------------------------------------------------

function [x1 d1] = alg2_Leg(Pn0,n) % Find the first root (note P_n'(0)==0)
% Approximate first root via asymptotic formula
k = ceil(n/2); theta = pi*(4*k-1)/(4*n+2);
x1 = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);

m = 30; % Number of terms in Taylor expansion.

% Recurrence Taylor coefficients (scaled & incl factorial terms).
M = 1/x1; % Scaling
zz = zeros(m,1); u = [Pn0 zeros(1,m)]; up = zeros(1,m+1); % Allocate storage
for k = 0:2:m-2
    up(k+2) = (k-n*(n+1)/(k+1))*u(k+1)/M^2;
    u(k+3) = up(k+2)/(k+2);
end
% Flip for more accuracy in inner product calculation.
u = u(m+1:-1:1); up = up(m+1:-1:1);

% % Note, terms are flipped for more accuracy in inner product calculation.
% zz = zeros(m,1); u = [zeros(1,m) Pn0]; up = zeros(1,m+1); % Allocate storage
% for k = 0:2:m-2
%     up(m-k) = (k-n*(n+1)/(k+1))*u(m-k+1)/M^2;
%     u(m-(k+1)) = up(m-k)/(k+2);
% end

% Newton iteration
x1k = ones(m+1,1); step = inf; l = 0;
while (abs(step) > eps) && (l < 10)
    l = l + 1;
    step = (u*x1k)/(up*x1k)/M;
    x1 = x1 - step;
    x1k = [1;cumprod(M*x1+zz)]; % Powers of h (This is the fastest way!)
    x1k = x1k(end:-1:1);        % Flip for more accuracy in inner product
end

% Get the derivative at this root, i.e. P'(x1).
d1 = M*(up*x1k);


