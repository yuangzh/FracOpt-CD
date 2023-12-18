function test
load aaa;
HandleObj = @(t)ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0);
[x1] = fminsearch(HandleObj,0); 
% if(HandleObj(x_best2)<HandleObj(x_best))
%     save aaa a2 a1 a0 b4 b3 b2 b1 b0
%     error('Not Best');
% end

x2 = nonconvex_prox_frac_l2l4(a2,a1,a0,b4,b3,b2,b1,b0,-inf,inf)

HandleObj(x1)-HandleObj(x2)

function [fobj,grad] = ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0)
up =  a2*t*t + a1*t + a0;
ooo = b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0;
if(ooo<0)
    fobj = 1e100;
    grad = 1e100;
    return;
end
down = sqrt(ooo);
fobj = up./down;
grad_up = 2*a2*t + a1;
grad_down = 0.5*(b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0)^(-0.5)*(4*b4*t*t*t + 3*b3*t*t + 2*b2*t + b1);
grad = grad_up - fobj.*grad_down;
grad = grad ./ down;