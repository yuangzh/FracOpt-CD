function [x,his,ts] = SpaseRecovery_ParametricCD(G,y,x0,lambda,k,speak,accuracy,MaxTime)
% min_x (0.5 ||Gx-y||_2^2  + lambda ||x||_1)  /  (lambda ||x||_{top_k})
accuracy = 1e-10;
initt = clock;

HandleObj = @(x)(0.5* norm(G*x-y)^2 + lambda * norm(x,1)) / (lambda*topksum(x,k));


x = x0;

n = length(x0);


gg = sum(G.*G,1)' + 1e-5;
Gy = G'*y;
Gx = G*x;

his = []; ts =[];
fobj = HandleObj(x);
his = [his;fobj];
ts = [ts;etime(clock,initt)];

last_k = 100;
stop_accuracy = accuracy;
changes = ones(last_k,1)*inf;
fobj_old = fobj;


for iter = 1:50000000
    %     i = randperm(n,1);
    
    
    [down,subgrad] = topksum(x,k);
    
    down = down*lambda;
    %     subgrad = subgrad * lambda;
    up = 0.5* norm(Gx-y)^2 +  lambda*norm(x,1) ;
    
    fobj = up / down;
    
    %     grad11 = G'*(G*x-y) - fobj*subgrad;
    %     ddd = abs(prox_l1(x-grad11/L,lambda/L)-x);
    %     [~,i] = min(ddd);
    
    
    i = mod(iter,n)+1;
%     ei = zeros(n,1);
%     ei(i) = 1;
    
    
    % 0.5 L t^2 + g_i t + lambda ||x+tei||_1 - c lambda ||x+tei||_{top_k}
    
    %     grad_i = G'*(Gx - y);
    %
    %     grad_i = grad_i(i);
    d = G(:,i);
    grad_i = Gx'*d - Gy(i);
%     HaneleObj = @(t) 0.5*gg(i)*t*t + grad_i*t + lambda*norm(x+t*ei,1)-fobj*lambda*topksum(x+t*ei,k);
    
    t = nonconvex_prox_diff_l1topk(gg(i),grad_i,x,i,k,lambda,fobj*lambda);
    % t = fminsearch(HaneleObj,0);
    x(i) = x(i)+t;
    Gx = Gx + t*G(:,i);
 
    if(speak)
    fprintf('iter:%d, fobj:%f\n',iter,fobj);
    end
 
    
    his(iter) = fobj;
    ts(iter) = etime(clock,initt);
    
 
    if etime(clock,initt) > MaxTime,
        break;
    end
    
    rel_change = abs((fobj - fobj_old)/max(1,fobj_old));
    changes = [changes(2:end);rel_change];
    fobj_old = fobj;
    
 
    if(mean(changes)<stop_accuracy),break;end
    
    
    
end

% plot(his)
% ddd
