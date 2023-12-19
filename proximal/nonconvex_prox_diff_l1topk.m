function best_t = nonconvex_prox_diff_l1topk(a,b,x,i,k,lambda,lambda2)
% min_t 0.5 a t^2 + b t + lambda |x_i+t| - lambda2 || x+tei ||_{top-k}

% x_i is in the top-k subset
% min_t 0.5 a t^2 + b t + lambda |x_i+t| - lambda2 |x+tei|
% grad = at + b + lambda*sign(x_i+t) - lambda2 sign(x_i+t) = 0
% (1) grad = at + b + lambda*(+1) - lambda2 (+1) = 0
t1 = (lambda2 - lambda - b)/a;
% (2) grad = at + b + lambda*(-1) - lambda2 (-1) = 0
t2 = (lambda - lambda2 - b)/a;

% x_i is in NOT the top-k subset
% min_t 0.5 a t^2 + b t + lambda |x_i+t|
% grad = at + b + lambda*sign(x_i+t) = 0
% (3) grad = at + b + lambda*(+1) = 0
t3 = (-lambda - b)/a;
% (4) grad = at + b + lambda*(-1) = 0
t4 = (lambda - b)/a;

t5 = -x(i);
ts = [t1;t2;t3;t4;t5]; 

for iii=1:length(ts)
    fs(iii) = ComputeObj(ts(iii),a,b,x,i,k,lambda,lambda2);
end


[~,ind]=min(fs);
best_t = ts(ind);

 
function [f] = ComputeObj(t,a,b,x,i,k,lambda,lambda2)
ei = 0*x;
ei(i) = 1;
x1 = x + t*ei;
[val]=sort(abs(x1),'descend');
f = 0.5*a*t^2 + b*t  + lambda*norm(x1,1) - lambda2*sum(val(1:k));


