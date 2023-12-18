function [c0,c1,c2,c3,c4] = ComputeCofPoly4(x,A,i,Ax)
% min_t f(t) ||A(x+tei)||_4^4
% c0 + c1 * t + c2*t*t + c3*t*t*t + c4*t*t*t*t
% norm(A*(x+t*ei),4)^4
% a = A*ei;

a = A(:,i);
b = Ax;
bb = b.*b;
bbbb = bb.*bb;
ab = a.*b;
aa = a.*a;
aaaa = aa.*aa;
aaab = ab.*aa;
abab = ab.*ab;
abbb = ab.*bb;
c0 = sum(bbbb);
c1 = sum(4*abbb);
c2 = sum(6*abab);
c3 = sum(4*aaab);
c4 = sum(aaaa);