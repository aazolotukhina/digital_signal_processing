function [x] = tikhonovSVD(U,s,V,b,lambda)
p=size(s,1);
f=s.^2./(s.^2 + lambda^2);
x = V(:,1:p)*(f.*(U(:,1:p)'*b)./s);
end