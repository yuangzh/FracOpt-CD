function best_t = nonconvex_prox_l1topk2(a,b,lambda1,lambda2,x,i,k)
% min_t 0.5 a t^2 + b t + lambda1 ||x+tei||_1 - lambda2 || x+tei ||_{top-k}

% case 1: x(i)+ t ÊôÓÚtop-k set
% at + b + lambda1|x_i+t| -  lambda2|x_i+t|  = 0 
% x_i + t = 0   =>   t = -x_i
% x_i + t > 0   =>   at + b + lambda1 -  lambda2 = 0   =>   t = (-b - lambda1 + lambda2) / a
% x_i + t < 0   =>   at + b - lambda1 +  lambda2 = 0   =>   t = (-b + lambda1 - lambda2) / a


% case 2: x(i)+ t ²»ÊôÓÚtop-k set
% at + b + lambda1 sign(x(i)+t) = 0
% x(i) + t = 0   =>   t = -x(i)
% x(i) + t > 0   =>   at + b + lambda1 = 0   =>   t = (-b - lambda1) / a
% x(i) + t < 0   =>   at + b - lambda1 = 0   =>   t = (-b + lambda1) / a

ts = [-x(i);(-b-lambda1 + lambda2)/a; (-b + lambda1 - lambda2)/a; (-b-lambda1)/a;(-b+lambda1)/a];
HandleObj = @(t)ComputeObj(t,a,b,lambda1,lambda2,x,i,k);
best_fobj = 1e100;
for i = 1:length(ts)
    fobjs(i) = HandleObj(ts(i));
end

[~,ind]=min(fobjs);
best_t = ts(ind);


 
function [f] = ComputeObj(t,a,b,lambda1,lambda2,x,i,k)
ei = 0*x;
ei(i) = 1;
r = x+t*ei;
[val,ind]=sort(abs(r),'descend');
% B = maxk(abs(r),k);
% sum(val(1:k)) - B
f = 0.5*a*t^2 + b*t + lambda1*abs(x(i)+t) - lambda2*sum(val(1:k));




