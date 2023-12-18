function [fobj,grad] = ComputeObjGradNorm4_sqrt(x,A)
% fobj = sqrt(||Ax||_4^4)
[fobj1,grad1] = ComputeObjGradNorm4(x,A);
fobj = sqrt(fobj1);
grad = 0.5 / sqrt(fobj1)*grad1;