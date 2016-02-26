function [z,p] = sphhankel_laplace(order)

B = cell(order+1,1);
z = cell(order+1,1);
p = cell(order+1,1);
for n=0:order
  % Recursion formula for nominator
  B{n+1} = zeros(1,n+1);
  for k=0:n-1
      B{n+1}(k+1) = ((2*n-k-1)*(2*n-k)) / (2*(n-k)) * B{n}(k+1);
  end
  B{n+1}(n+1) = 1;
end

for n=0:order
  % Flip nominator polynoms
  B{n+1} = B{n+1}(end:-1:1);
  % Zeros
  z{n+1} = roots(B{n+1});
  % Poles (are always zero)
  p{n+1} = zeros(order,1);
end

z = z{order+1};
p = p{order+1};
