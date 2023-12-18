function test_nonconvex_prox_frac_l2l4
clc;clear all;close all;

% min A / B

% (a2 t^2 + a1 t + a0)^2
% --------------------------------------
% (b4 t^4 + b3 t^3 + b2 t^2 + b1 t^1 + b0)
for iter = 989:10000000000000
rand('seed',iter);
randn('seed',iter);

a2 = rand(1); a1 = randn(1); a0 = rand(1);
%[a0,a1,a2] =  GenMeaningfulCof0;
[b0,b1,b2,b3,b4] = GenMeaningfulCof;

a2 = 1;
a1 = 2.0000;
a0 = 1;
b4 =  1.1803e+13;
b3 =   4.7213e+13;
b2 =   7.0819e+13;
b1 = 4.7213e+13;
b0 = 1.1803e+13;



HandleObj = @(t)ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0);

L = -1000;
U = 1000;
% A' B = A B'
% 2 (a2 t^2 + a1 t + a0) (2 a2 t + a1) o (b4 t^4 + b3 t^3 + b2 t^2 + b1 t^1 + b0)
% =
% (a2 t^2 + a1 t + a0)^2 o ( 4 b4 t^3 + 3 b3 t^2 + 2 b2 t + b1)

% 2 (2 a2 t + a1) o (b4 t^4 + b3 t^3 + b2 t^2 + b1 t^1 + b0)
% =  (a2 t^2 + a1 t + a0) o ( 4 b4 t^3 + 3 b3 t^2 + 2 b2 t + b1)

 
 
 








[x1] = solve_by_matlab(a2,a1,a0,b4,b3,b2,b1,b0,L,U);
[x2] = nonconvex_prox_frac_l2l4(a2,a1,a0,b4,b3,b2,b1,b0,L,U);
[f1,g1] = HandleObj(x1);
[f2,g2] = HandleObj(x2);
fprintf('matlab: %f %f %f\n',x1,f1,g1);
fprintf('  ours: %f %f %f\n',x2,f2,g2);
fprintf('  diff: %d, %e\n',iter,f2-f1);

HandleObj(0)
HandleObj(1)
 
 ddd
t = x1;
up =  a2*t*t + a1*t + a0;
down = sqrt(b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0);

 
 
if (f2 - f1>  1e-2 * max(abs([1 f1 f2])))
    ddd
end
 
end

function [x] = solve_by_matlab(a2,a1,a0,b4,b3,b2,b1,b0,L,U)
HandleObj = @(t)ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0);
x = fminbnd(HandleObj,L,U);
% x = fminsearch(HandleObj,x);

% xs = [-2:0.01:-0.5];
% for i=1:length(xs)
%     [fs(i),gs(i)] = HandleObj(xs(i));
% end
%  
% plot(xs,fs); 
% figure;plot(xs,gs)
% % t =1.234;

function [a0,a1,a2] = GenMeaningfulCof0
% || x + t ei||_2^2

n = 1000;
x = randn(n,1);
a0 = x'*x;
i = randperm(n,1);
a1 = 2*x(i);
a2 = 1;

function  [c0,c1,c2,c3,c4] = GenMeaningfulCof
mn = randperm(1000,2);
m = mn(1);
n = mn(2);


A = randn(m,n);
x = randn(n,1);
i = randperm(n,1);
[c0,c1,c2,c3,c4] = ComputeCofPoly4(x,A,i);



function [fobj,grad] = ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0)

ooo = b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0;
if(ooo<0)
    fobj = 1e100;
    grad = 1e100;
    return;
end

up =  a2*t*t + a1*t + a0;
down = sqrt(ooo);
fobj = up./down;
grad_up = 2*a2*t + a1;
grad_down = 0.5*(b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0)^(-0.5)*(4*b4*t*t*t + 3*b3*t*t + 2*b2*t + b1);
grad = grad_up - fobj.*grad_down;
grad = grad ./ down;

% 2*a2*t + a1 = (a2*t*t + a1*t + a0) / sqrt(b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0)  o  0.5*(b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0)^(-0.5)*(4*b4*t*t*t + 3*b3*t*t + 2*b2*t + b1)
% 2*(2*a2*t + a1) (b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0) = (a2*t*t + a1*t + a0) (4*b4*t*t*t + 3*b3*t*t + 2*b2*t + b1)
% 2*(2*a2*t + a1) (b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0) = (a2*t*t + a1*t + a0) (4*b4*t*t*t + 3*b3*t*t + 2*b2*t + b1)
