function test_inequality
clc;
while 1
m = 11;
n = 11;
A = randn(m,n);
b = randn(m,1);

x = randn(n,1);
t = randn(1);
i = randperm(n,1);
ei = zeros(n,1);
ei(i) = 1;

c = diag(A'*A);

[fx,gx] = ComputeObj(x,A,b);
y = x + t * ei;
[fy,gy] = ComputeObj(y,A,b);

tt1 = fy 
tt2 = fx + t*gx(i) + c(i)*t*t/2

if(tt1-tt2>1e-8)
    tt1 - tt2
    
    ddd
end

end

function [fobj,grad] = ComputeObj(x,A,b)
fobj = 0.5*norm(max(0,A*x-b))^2;
grad = A'*max(0,A*x-b);
