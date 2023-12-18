function test_nonconvex_prox_frac_l1
addpath('util');
for iter = 1:1000000,
seed = iter
rand('seed',seed);
randn('seed',seed);
n = 31;
a = randn(1);
b = randn(1);
c = randn(1);
g = randn(n,1);
d = randn(n,1);


HandleObj = @(t) (0.5*a*t*t + b*t + c) / (norm(g*t+d,1)^2);
t1 = fminsearch(HandleObj,0);
t2 = nonconvex_prox_frac_l1(a,b,c,g,d);


if(HandleObj(t1) - HandleObj(t2)<-0.0001)
    iter
    dd
end

end
