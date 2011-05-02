function beta_nk=beta1(n,k)
if k==n
    beta_nk=1;
else
    beta_nk=(((2.*n-k-1).*(2.*n-k) )./(2.*(n-k))).*beta1(n-1,k);
end