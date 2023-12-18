function [x_best,his]=AccerlatedProximalGradient(x,HandleObjSmooth,HandleObjNonSmooth,HandleProx,maxiter,L)
% This program solves the following optimization problem:
% f(x) + g(x)
% where we assume that f is smooth g is non-smooth
% HandleObjSmooth:           x   ->  [fobj,grad]
% HandleObjNonSmooth:        x   ->  [fobj]
% HandleProx:          [theta,a] ->  arg min_{x} 0.5 theta || x - a ||^2 + g(x)

% One xample:

% function example_LeastR
% clear, clc;
% %  min  1/2 || A x - y||^2 + lambda * ||x||_1
% % f(x) + g(x)
% % 0.5 L ||x-xt||^2 + <x-xt,g> + g(x)
% % 0.5 L ||x-(xt-g/L)||^2 + g(x)
% % proximal mapping:
% 
% % 0.5 theta ||x - a||^2 + g(x)
% 
% 
% m=1000;  n=100;    % The data matrix is of size m x n
% A=randn(m,n);       % the data matrix
% y = randn(m,1);
% lambda=0.2;
% HandleObjSmooth = @(x)computeObj(x,A,y);
% HandleObjNonSmooth = @(x)lambda*sum(abs(x));
% x=zeros(n,1);
% HandleProx = @(theta,a)computeprox(theta,a,lambda);
% [x1, his]= AccerlatedProximalGradient(x,HandleObjSmooth,HandleObjNonSmooth,HandleProx);
% plot(his)
% 
% function [fobj,grad] = computeObj(x,A,y)
% diff = A*x-y;
% fobj = 1/2*norm(diff)^2 ;
% grad = A'*diff ;
% 
% function [x] = computeprox(theta,a,lambda)
% % 0.5 theta ||x - a||^2 + g(x)
% [x] = threadholding_l1(a,lambda/theta);
% 
% function [x] = threadholding_l1(a,lambda)
% % solving the following OP:
% % min_{x} 0.5 ||x - a||^2 + lambda * sum(abs(x))
% x = sign(a).*max(0,abs(a)-lambda);

% last modified: 2016-01-29

[n,d]=size(x);
flag=0;

xp=x;
xxp=zeros(n,d);
alpha=1; s=x;
x_best = x;
f_best = HandleObjSmooth(x)+ HandleObjNonSmooth(x);
for iterStep=1:maxiter,
%     fprintf('*');
    [f_old,g_old] = HandleObjSmooth(s);
    xp=x;
    max_in = 10000;
    for in=1:max_in,
        v=s-g_old/L;
        % min_{x} 0.5 L ||x-(xt-g/L)||^2 + z * ||x||_1
        [x]=HandleProx(L,v);
        v=x-s;
        r_sum=mdot(v,v);
        f_new = HandleObjSmooth(x);
        l_sum = f_new - f_old - mdot(g_old,x-s);
        if (r_sum <=1e-100)
            flag=1; % this shows that, the gradient step makes little improvement
            break;
        end
        if(l_sum <= 0.5*r_sum*L)
            break;
        else
%             L=max(2*L,l_sum/r_sum);
            L=2*L;
        end
        if(in==max_in)
            fprintf('warning! Lipschitz too large!');
        end
    end
    alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;
    beta=(alphap-1)/alpha;    s=x + beta* xxp;
    xxp=x-xp;
    f_curr = f_new  + HandleObjNonSmooth(x);
    if(f_curr<f_best)
        f_best = f_curr;
        x_best = x;
    end
    his(iterStep) = f_curr;
%     if (iterStep>50 && flag),break;end
end
his = his(:);

function [r] = mdot(x,y)
r = sum(sum(x.*y));