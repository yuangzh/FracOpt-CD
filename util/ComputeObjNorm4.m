function [fobj] = ComputeObjNorm4(x,A)
x = norm_one(x);
fobj = norm(A*x,4);