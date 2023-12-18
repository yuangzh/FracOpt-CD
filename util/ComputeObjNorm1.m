function [fobj] = ComputeObjNorm1(x,A)
% x = norm_one(x);
% fobj = norm(A*x,1);
fobj = norm(x)  / norm(A*x,1);
