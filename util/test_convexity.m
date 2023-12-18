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
    t = rand(1);
    z = t * x + (1-t)*y;


    F = @(x,num)ComputeObj(x,A,b,n,num);
    f = @(x,num)ComputeObj_up(x,A,b,n,num);
    g = @(x,num)ComputeObj_down(x,A,b,n,num);

%  
%     ooo = f(z,1);
%     ooo = t*f(x,1) + (1-t)*f(y,1);
    

    
%     tttt = ooo/ (t*g(x,1) + (1-t)*g(y,1)) - max(F(x,1),F(y,1));
   tttt = F(z,1) - (t*F(x,1)+(1-t)*F(y,1));
 
     
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