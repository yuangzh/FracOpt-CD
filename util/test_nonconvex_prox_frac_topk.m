function test_nonconvex_prox_frac_topk
clc;clear all;close all;
% min_t (0.5 a t^2 + b t + c) / || x+tei ||_{top-k}

for iter = 1:100000000
    iter
    rand('seed',iter);
    randn('seed',iter);
    m = 115;
    a = rand(1)*1000*rand(1);
    b = randn(1)*1000*rand(1);
    c = randn(1)*1000*rand(1);
    x = randn(m,1)*1000*rand(1);
    i = randperm(m,1);
    k = randperm(m,1);


    b = b*max(0,randn(1));
    x = x.*max(0,randn(m,1));
    

    HandleObj = @(t)ComputeObj(t,a,b,c,x,i,k);
    x1 = fminsearch(HandleObj,0);

    x2 = nonconvex_prox_frac_topk(a,b,c,x,i,k);

    f1 = HandleObj(x1);
    f2 = HandleObj(x2);
    fprintf('iter:%f, fobj:%.5e %.5e, f1-f2:%f\n',iter,f1,f2,f1-f2);
 
    if(f2>f1 + 1e-5*abs(mean([f1;f2])))
        ei = zeros(m,1);
        ei(i)  = 1;
     
        i
        k
        x+x1*ei
        f1
        f2
        x1
        x2
        t = x1;
 
        dddd
    end
    
    
end

 

 
function [f] = ComputeObj(t,a,b,c,x,i,k)
ei = 0*x;
ei(i) = 1;
[val]=sort(abs(x+t*ei),'descend');
% B = maxk(abs(r),k);
% sum(val(1:k)) - B
f = (0.5*a*t^2 + b*t + c) /  sum(val(1:k));

