function B = sphexpS_mono_multiscatter(A, xq, R, sigma, f, conf)
  
%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargmatrix(A,xq);
isargpositivescalar(f);
isargvector(R, sigma);
if nargin<nargmax
  conf = SFS_config;
else
  isargstruct(conf);
end
if length(R) == 1
  R = repmat(R,[1, size(A,2)]);
elseif length(R) ~= size(A,2)
  error('%s: Length of R does not match size of A',upper(mfilename));
end
if length(sigma) == 1
  sigma = repmat(sigma,[1, size(A,2)]);
elseif length(sigma) ~= size(A,2)
  error('%s: Length of R does not match size of A',upper(mfilename));
end

%% ===== Configuration ==================================================
Nse = conf.scattering.Nse;

%% ===== Computation ====================================================
if size(A,1) ~= (Nse + 1)^2
  error('%s: size of A does not match Nse',upper(mfilename));
end

Nl = size(A,1);
Nq = size(A,2);
L = zeros(Nq*Nl);

E = ones(Nl,1);
for qdx=1:Nq
  selectq = ((qdx-1)*Nl+1):(qdx*Nl);
  L(selectq,selectq) = diag(1./sphexpS_mono_scatter(E, R(qdx), sigma(qdx), f, conf));  
end

for qdx=1:Nq
  selectq = ((qdx-1)*Nl+1):(qdx*Nl);
  for pdx=(qdx+1):Nq
    selectp = ((pdx-1)*Nl+1):(pdx*Nl);    
    [SRpq, SRqp] = sphexpSR_mono(xq(qdx,:)-xq(pdx,:), f, conf);
    L(selectp,selectq) = -SRpq;
    L(selectq,selectp) = -SRqp;
  end
end

A = reshape(A,[],1);

B = L\A;
B = reshape(B,Nl,Nq);
