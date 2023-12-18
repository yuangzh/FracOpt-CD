function test_nonconvex_diff_frac_l2l4

% min_{L<=t<=U} (a2 t^2 + a1 t + a0)^2 - C £¨b4 t^4 + b3 t^3 + b2 t^2 + b1 t^1 + b0£©

clc;clear all;close all;

% min A / B

for iter = 989:10000000000000
rand('seed',iter);
randn('seed',iter);

% a2 = rand(1); a1 = randn(1); a0 = rand(1);
[a0,a1,a2] =  GenMeaningfulCof0;
[b0,b1,b2,b3,b4] = GenMeaningfulCof;
C = rand(1);  
L = -1000;
U = 1000;
t = randn(1);
ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0,C)
ddddd
HandleObj = @(t)ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0,C);

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







function x_best = nonconvex_prox_diff_l2l4(a2,a1,a0,b4,b3,b2,b1,b0,L,U)
% min_{L<=t<=U} a2 t^2 + a1 t + a0 - C sqrt£¨b4 t^4 + b3 t^3 + b2 t^2 + b1 t^1 + b0£©

% We assume that the denominator is positive
% a2,a1,a0,b4,b3,b2,b1,b0,L,U

Proj = @(x) min(U,max(L,x));
HandleObj = @(t)ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0);
cof_5 = 4*a2.*b4;
cof_4 = 4*a2.*b3 + 2*a1.* b4;
cof_3 = 4*a2.*b2 + 2*a1.*b3; 
cof_2 = 4*a2.*b1+ (2*a1) .* ( b2);
cof_1 = 4*a2.*b0  + (2*a1) .* (b1);
cof_0 = 2*a1.*b0;
% left = cof_5 * t^5 + cof_4*t^4 + cof_3*t^3 + cof_2*t^2 +  cof_1*t   + cof_0;




cof2_5 = (a2) .* ( 4* b4 );
cof2_4 = (a2) .* (3* b3) +    ( a1 ) .* ( 4* b4);
cof2_3 = (a2) .* (2 *b2 ) +    ( a1) .* ( 3* b3 )  +   ( a0) .* ( 4* b4 ); 
cof2_2 = (a2 ) .*  b1 +    ( a1) .* ( 2 *b2 )  +   ( a0) .* ( 3* b3);
cof2_1 =   ( a1) .* ( b1)  +   ( a0) .* ( 2 *b2 );
cof2_0 =   ( a0) .* ( b1);

Left = @(t)(4.*a2.*t + 2.*a1) .* (b4 *t.^4 + b3.* t.^3 + b2.*t.^2 + b1.*t + b0);
Right = @(t)(a2.*t.^2 + a1.*t + a0) .* ( 4.* b4.* t.^3 + 3.*b3.*t.^2 + 2 *b2.*t + b1);

c4 = cof_4-cof2_4; c3 = cof_3-cof2_3; c2 = cof_2-cof2_2; c1 = cof_1-cof2_1; c0 = cof_0-cof2_0;


% find t such that left == right
%  cof_5 * t^5 + cof_4*t^4 + cof_3*t^3 + cof_2*t^2 +  cof_1*t   + cof_0  = cof2_5 * t^5 + cof2_4*t^4 + cof2_3*t^3 + cof2_2*t^2 +  cof2_1*t   + cof2_0
% it reduces to 
%  cof_4*t^4 + cof_3*t^3 + cof_2*t^2 +  cof_1*t   + cof_0  = cof2_4*t^4 + cof2_3*t^3 + cof2_2*t^2 +  cof2_1*t   + cof2_0

xs1 = roots([c4;c3;c2;c1;c0]);
xs2 = roots([b4;b3;b2;b1;b0]);

xs1 = [abs(xs1);real(xs1)];
xs2 = [abs(xs2);real(xs2)];

xs = [0;xs1;xs2];


xs = [L;U; Proj(xs)];

for iii=1:length(xs)
fs(iii) = HandleObj(xs(iii));
end

[~,ind] = min(fs);
x_best = xs(ind);

% ind = find(fs==fs(ind));
% x_candidate = xs(ind);
% [~,ind2] = min(abs(x_candidate));
% x_best = x_candidate(ind2);





function [fobj] = ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0,C)
O = b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0;

fobj = (a2*t*t + a1*t + a0)^2 - C*O;


cof_4 = a2*a2;
cof_3 = 2*a2*a1;
cof_2 = 2*a2*a0 + a1*a1;
cof_1 = 2*a1*a0;
cof_0 = a0*a0;

cof_4 = cof_4 

