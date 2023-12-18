function [x] = prox_l2(a,lambda)
% solving the following OP:
% min_{x} 0.5||x-a||^2 + lambda * ||x||_2
% max_y min_x 0.5||x-a||^2 + lambda <x,y>, s.t. ||y||<=1
% grad_x = x - a + lambda y = 0   =>    x = a - lambda y
% max_y 0.5||a - lambda y -a||^2 + lambda <a - lambda y,y>, s.t. ||y||<=1
% max_y 0.5lambda|| y||^2 +  <a - lambda y,y>, s.t. ||y||<=1
% max_y 0.5lambda|| y||^2 + <a,y> - lambda <y,y>, s.t. ||y||<=1
% max_y -0.5lambda|| y||^2 + <a,y>, s.t. ||y||<=1
% min_y 0.5lambda||y||^2 - <a,y>, s.t. ||y||<=1
% min_y 0.5||y||^2 - <a/lambda,y>, s.t. ||y||<=1
% min_y 0.5||y - a/lambda||^2, s.t. ||y||<=1
y = projectL2(a/lambda,1);
x = a - lambda*y;


function [w] = projectL2(w,tau)
nw = norm(w,'fro');
if nw > tau
    w = w*(tau/nw);
end


