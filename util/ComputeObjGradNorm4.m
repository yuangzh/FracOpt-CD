function [fobj,grad] = ComputeObjGradNorm4(x,A)
% fobj = ||Ax||_4^4
[m,n] = size(A);
Ax = A*x;
Ax3 = Ax.*Ax.*Ax;
fobj = norm(Ax,4)^4;
% grad = zeros(n,1);
% for i=1:m,
% a = A(i,:)';
% one_grad = 4*Ax(i)^3*a;
% grad = grad + one_grad;
% end
grad = 4*sum(diag(Ax3)*A,1)';