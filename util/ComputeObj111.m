function [cof1,cof2,cof3,cof4] = ComputeObj111(x,A,c,i)
% min_t f(t) norm(x+tei)^4 - c ||A(x+tei)||_4^4  + 0.5*beta*t^2;
% a = A*ei;
a = A(:,i);
b = A*x;
bb = b.*b;
bbbb = bb.*bb;
ab = a.*b;
aa = a.*a;
aaaa = aa.*aa;
aaab = ab.*aa;
abab = ab.*ab;
abbb = ab.*bb;
xx = x'*x;

cof_0 = xx*xx;
cof_1 = 4*xx*x(i);
cof_2 = 2*xx+4*x(i)^2;
cof_3 = 4*x(i);
cof_4 = 1;

cof_00 = bbbb;
cof_11 = 4*abbb;
cof_22 = 6*abab;
cof_33 = 4*aaab;
cof_44 = aaaa;

cof0 = cof_0  - c*sum(cof_00);
cof1 = cof_1  - c*sum(cof_11);
cof2 = cof_2  - c*sum(cof_22);
cof3 = cof_3  - c*sum(cof_33);
cof4 = cof_4  - c*sum(cof_44);

% fobj =   norm(x+ei*t)^4 - c*norm(t*a+b,4)^4 ;
% fobj2 =  cof0 + cof1*t + cof2*t*t + cof3*t*t*t + cof4*t*t*t*t;

% fobj2 - fobj