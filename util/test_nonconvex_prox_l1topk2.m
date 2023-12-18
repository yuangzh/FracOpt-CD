function test_nonconvex_prox_l1topk2
clc;clear all;close all;
% min_t 0.5 a t^2 + b t + lambda1 ||x+tei||_1 - lambda2 || x+tei ||_{top-k}



for iter = 15:100000000
    iter
    rand('seed',iter);
    randn('seed',iter);
    m = 10;
    a = rand(1)*100*rand(1);
    b = randn(1)*100*rand(1);
    x = randn(m,1)*100*rand(1);
    lambda1 = rand(1)*100*rand(1);
    lambda2 = rand(1)*100*rand(1);
    i = randperm(m,1);
    k = randperm(m,1);
    

    b = b*max(0,randn(1));
    x = x.*max(0,randn(m,1));
    lambda1 = lambda1.*max(0,randn(1));
    lambda2 = lambda2.*max(0,randn(1));
    

    HandleObj = @(t)ComputeObj(t,a,b,lambda1,lambda2,x,i,k);
    x1 = fminsearch(HandleObj,0);
 
    x2 = nonconvex_prox_l1topk2(a,b,lambda1,lambda2,x,i,k);

    f1 = HandleObj(x1);
    f2 = HandleObj(x2);
    fprintf('iter:%f, fobj:%.5e %.5e\n',iter,f1,f2);
    
    if(f2>f1 + 1e-6*abs(mean([f1;f2])))
        f1
        f2
        x1
        x2
        dddd
    end
    
    
end

 

 
function [f] = ComputeObj(t,a,b,lambda1,lambda2,x,i,k)
ei = 0*x;
ei(i) = 1;
r = x+t*ei;
[val,ind]=sort(abs(r),'descend');
% B = maxk(abs(r),k);
% sum(val(1:k)) - B
f = 0.5*a*t^2 + b*t + lambda1*abs(x(i)+t) - lambda2*sum(val(1:k));




