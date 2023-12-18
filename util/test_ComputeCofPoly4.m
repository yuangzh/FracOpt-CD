clc; clear all; close all;

m = 20;
n = 11;
A = randn(m,n);
i = randperm(n,1);
ei = zeros(n,1);
ei(i) = 1;

x = randn(n,1);
alpha = 1111.2345;


[b0,b1,b2,b3,b4] = ComputeCofPoly4(x,A,i);
sqrt(b4*alpha^4 + b3*alpha^3 + b2*alpha^2 + b1*alpha^1 + b0) - norm(A*(x+alpha*ei),4)^2

 

