function test_nonconvex_prox_frac_linf
% min_t (at^2 + bt + c) / ||gt+d||_inf^2
% We assume that (at^2 + bt + c)>=0 for all t, otherwise 
% this problem cannot be solved exactly usning breakpoint search method


addpath('util');
for iter = 1:1000000,
    seed = iter
    rand('seed',seed);
    randn('seed',seed);
    n = 2;
    a1 = randn(1);
    b1 = randn(1);
    
    g = randn(n,1);
    d = randn(n,1);
    
    
%     a = a1*a1;
%     b = 2*a1*b1;
%     c = b1*b1;
    
    a = randn(1);
    b = randn(1);
    c = randn(1);
    
    % (at+b)^2 = aatt + 2abt + bb
    
    HandleObj = @(t) (a*t*t + b*t + c) / (norm(g*t+d,inf)^2);
    t1 = fminsearch(HandleObj,0);
    t2 = nonconvex_prox_frac_linf(a,b,c,g,d);
    
    
    if(HandleObj(t1) - HandleObj(t2)<-0.0001)
        
        HandleObj(t1)
        HandleObj(t2)
        p = [g;-g]; w = [d;-d]; 
        [fobj,grad] = ComputeObj_t(t1,a,b,c,p,w);
        fobj
        grad
        %
        %     t1
        %     t2
        %
        ts = [-5:0.0001:5];
        fffffs = [];
        for ii=1:length(ts)
            fffffs(ii) = HandleObj(ts(ii));
        end
        plot(ts,fffffs)
        
        
        dddddddddd
    end
    
end
function [fobj,grad] = ComputeObj_t(t,a,b,c,p,w)
% (a t t + b t + c ) / (max(pt+w)^2)
o = p*t+w;
[val,ind] = max(o);
up = a*t*t + b*t + c;
down = (val^2);
grad_up = 2*a*t + b;
grad_down = 2*val*p(ind);
 
fobj = up / down;
grad = (grad_up * down - up*grad_down) / (down^2);
 up
 down
 o
% grad2 = grad_up / down - up*grad_down / (down^2);

 
% grad =  ( 2*a*w(ind) - b*2*p(ind) + b * p(ind) ) *t +        b * w(ind)  - c*2*p(ind)


