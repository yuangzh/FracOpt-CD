function [x,his,ts] = SpaseRecovery_FractionalCD(G,y,x0,lambda,k,MaxTime)
% min_x (0.5 ||Gx-y||_2^2  + lambda ||x||_1)  /  ||x||_{top_k}



initt = clock;
x = x0;
his = []; ts = [];
n = length(x0);



gg = sum(G.*G,1)';
Gy = G'*y;
Gx = G*x;
const = 0.5*y'*y;
HandleSquare = @(x) 0.5*norm(G*x-y)^2 + lambda*norm(x,1);
fobj = HandleSquare(x) + topksum(x,k);


% his = [his;fobj];
ts = [ts;etime(clock,initt)];

HandleObj = @(x)HandleSquare(x) / topksum(x,k);

last_k = 100;
stop_accuracy = 1e-14;
changes = ones(last_k,1)*inf;
fobj_old = 1e10;
for iter = 1:4694
    grad = G'*(G*x-y);
    i = randperm(n,1);
 

    c = HandleObj(x);
%     i = mod(iter,n)+1;
    ei = zeros(n,1);
    ei(i) = 1;

    
% 0.5 ||Gx-y||_2^2  + lambda ||x||_1  - c ||x||_{top_k}
    
o = G*x - y;
w = G*ei;
% min_t (   0.5 ||G(x+tei)-y||_2^2  + lambda ||(x+tei)||_1)  /  ||(x+tei)||_{top_k}
% min_t (   0.5 || o + t w ||_2^2  + lambda ||(x+tei)||_1)  /  ||(x+tei)||_{top_k}
% min_t (   0.5 w'w t^t + <o,w> t + 0.5 o'o  + lambda ||(x+tei)||_1)  /  ||(x+tei)||_{top_k}
 
t = nonconvex_prox_frac_l1topk(w'*w,o'*w,0.5*o'*o,x,i,k,lambda);
 

    x(i) = x(i)+t;
    fprintf('iter:%d, fobj:%f\n',iter,c);
    his = [his;c];

end
plot(his)
ddddd