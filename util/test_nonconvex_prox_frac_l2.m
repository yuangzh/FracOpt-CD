function test_nonconvex_prox_frac_l2
clc;clear all;close all;
% min_t (0.5 a t^2 + b t + c) / || x+tei ||_2

for iter = 1:100000000
    iter
    rand('seed',iter);
    randn('seed',iter);
    m = 15;
    a = rand(1)*1000*rand(1);
    b = randn(1)*1000*rand(1);
    c = randn(1)*1000*rand(1);
    x = randn(m,1)*1000*rand(1);
    i = randperm(m,1);

    

    HandleObj = @(t)(0.5*a*t^2 + b*t + c) /  sqrt(x'*x + 2*x(i)*t + t*t);
    x1 = fminsearch(HandleObj,0);

    x2 = nonconvex_prox_frac_l2(a,b,c,x,i);

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

 

 






