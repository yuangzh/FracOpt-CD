function test_nonconvex_prox_diff_l1topk
clc;clear all;close all;
% min_t 0.5 a t^2 + b t + lambda ||x+tei||_1 - lambda2 || x+tei ||_{top-k}

for iter = 1:100000000
 
    rand('seed',iter);
    randn('seed',iter);
    m = 115;
    
    a = rand(1)*1000*rand(1);
    b = randn(1)*1000*rand(1);
    x = randn(m,1)*1000*rand(1);
    i = randperm(m,1);
    k = randperm(m,1);
    lambda = rand(1)*1000*rand(1);
    lambda2 = rand(1)*1000*rand(1);


    b = b*max(0,randn(1));
    x = x.*max(0,randn(m,1));
    

    HandleObj = @(t)ComputeObj(t,a,b,x,i,k,lambda,lambda2);
    x1 = fminsearch(HandleObj,0);

    x2 = nonconvex_prox_diff_l1topk(a,b,x,i,k,lambda,lambda2);
 
    f1 = HandleObj(x1);
    f2 = HandleObj(x2);
    fprintf('iter:%f, fobj:%.5e %.5e, f1-f2:%f\n',iter,f1,f2,f1-f2);
 
    if(f2>f1 + 1e-2*abs(mean([f1;f2])))
 
 
        dddd
    end
    
    
end

 

 
function [f] = ComputeObj(t,a,b,x,i,k,lambda,lambda2)
ei = 0*x;
ei(i) = 1;
x1 = x + t*ei;
[val]=sort(abs(x1),'descend');
% B = maxk(abs(r),k);
% sum(val(1:k)) - B
f = 0.5*a*t^2 + b*t  + lambda*norm(x1,1) - lambda2*sum(val(1:k));





































