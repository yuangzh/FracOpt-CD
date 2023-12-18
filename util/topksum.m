function [fobj,subgrad] = topksum(x,k)
% f(x) = top_k(abs(x))
n = length(x);
absx = abs(x);
[sort_absx,ind]=sort(absx,'descend');
fobj = sum(sort_absx(1:k));
subgrad = zeros(n,1);
subgrad(ind(1:k))=1;
subgrad = mysign(x).*subgrad;


function [r] = mysign(x)
r = sign(x);
I = find(x==0);
r(I) = sign(randn(length(I),1));