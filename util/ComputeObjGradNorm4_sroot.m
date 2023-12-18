function [fobj,grad] = ComputeObjGradNorm4_sroot(x,A)
% fobj = (||Ax||_4^4)^(0.5)

[m,n] = size(A);
Ax = A*x;
Ax3 = Ax.*Ax.*Ax;
fobj = sqrt(norm(Ax,4)^4);
% grad = zeros(n,1);
% for i=1:m,
% a = A(i,:)';
% one_grad = 4*Ax(i)^3*a;
% grad = grad + one_grad;
% end
grad = 0.5*(norm(Ax,4)^4)^(-0.5)*4*sum(diag(Ax3)*A,1)';