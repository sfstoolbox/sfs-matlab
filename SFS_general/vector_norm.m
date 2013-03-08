function y = vector_norm(x,dim)
y = sum(abs(x).^2,dim).^(1/2);