function [z,p] = sphhankel_laplace(order)

hankel = struct(...
  'zeros',[],...
  'poles',[],...
  'scale',[],...
  'delay',[],...
  'nominator',[],...
  'denominator',[]);

hankel = repmat(hankel,order+1,1);
for n=0:order
  % nominator
  hankel(n+1).nominator = zeros(1, n+1);
  for k=0:n-1
    hankel(n+1).nominator(k+1) = ...
      ((2*n-k-1)*(2*n - k))/(2*(n-k))*hankel(n).nominator(k+1);
  end
  hankel(n+1).nominator(n+1) = 1;
end

for n=0:order
  % flip nominator polynoms
  hankel(n+1).nominator = hankel(n+1).nominator(end:-1:1);
  % zeros
  hankel(n+1).zeros = roots(hankel(n+1).nominator);
  % denominator
  hankel(n+1).denominator = zeros(1, n+2);
  hankel(n+1).denominator(1) = 1;
  % poles
  hankel(n+1).poles = roots(hankel(n+1).denominator);
end

z = hankel(order).zeros;
p = hankel(order).poles;
