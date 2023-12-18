function test_convexity
clc;clear all;close all;
% f(x) / g(x)
% f(x) is convex
% g(x) is concave

m = 1;
n = 11;


while 1
    
    A = randn(m,n);
    b = randn(m,1);

    Proj = @(x) max(-1,min(x,1));
    x = 1*randn(n,1);
    y = 1*randn(n,1);
    x = Proj(x);
    y = Proj(y);

    
    
    F = @(x,num)ComputeObj(x,A,b,n,num);
    f = @(x,num)ComputeObj_up(x,A,b,n,num);
    g = @(x,num)ComputeObj_down(x,A,b,n,num);

 
    tttt = f(x,1)/g(x,1) + F(x,2)'*(y-x) - f(y,1)/g(y,1);
 
 
    tttt = f(x,1)/g(x,1) + ((f(x,2) - F(x,1)*g(x,2)) )'*(y-x)/ g(x,1) - f(y,1)/g(y,1);

    tttt = f(x,1)*g(y,1) + g(y,1)*(  ( f(x,2) - F(x,1)*g(x,2)) )'*(y-x) - f(y,1)*g(x,1);
 
    tttt2 = f(x,1)*g(y,1)   +     g(y,1)*f(x,2)'*(y-x)     -     g(y,1)*F(x,1)*g(x,2)'*(y-x)             - f(y,1)*g(x,1);
    tttt = f(x,1)*g(y,1)   +     g(y,1)*f(x,2)'*(y-x)     -     g(y,1)*F(x,1)*(g(y,1)-g(x,1))             - f(y,1)*g(x,1);
 
    tttt = f(x,1)*g(y,1)   +     g(y,1)*f(x,2)'*(y-x)     -     g(y,1)*f(x,1)/g(x,1)*(g(y,1)-g(x,1))             - f(y,1)*g(x,1);

    tttt =     g(y,1)*g(x,1)*f(x,2)'*(y-x)     +    f(x,1)*(2*g(x,1)*g(y,1)   -      g(y,1)*g(y,1)           )          - f(y,1)*g(x,1)*g(x,1);
     
      fprintf('   convex: %.2f\n',tttt)
    if(  tttt > 0)
        fprintf('nonconvex: %.2f\n',tttt)
    fgggsdd
    end
    
end;

function out = ComputeObj_up(x,A,b,n,num)
fobj = 0.5*norm(A*x-b)^2;
grad = A'*(A*x-b);
if(num==1)
    out = fobj;
else
    out = grad;
end


function out = ComputeObj_down(x,A,b,n,num)
fobj = (n - norm(x,1));
grad = - sign(x);
if(num==1)
    out = fobj;
else
    out = grad;
end




function out = ComputeObj(x,A,b,n,num)
up = 0.5*norm(A*x-b)^2;
down = (n - norm(x,1));
grad_up = A'*(A*x-b);
grad_down = - sign(x);
fobj = up / down;
grad = grad_up - fobj * grad_down;
grad = grad / down;

if(num==1)
    out = fobj;
else
    out = grad;
end