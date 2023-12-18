function [x,his,ts] = SpaseRecovery_DPA(G,y,x,lambda,k,speak,accuracy,MaxTime)
% min_x (0.5 ||Gx-y||_2^2  + lambda ||x||_1)  /  (lambda ||x||_{top_k})
initt = clock;
HandleObj = @(x)(0.5* norm(G*x-y)^2 + lambda * norm(x,1)) / (lambda*topksum(x,k));


his = []; ts =[];
fobj = HandleObj(x);
his = [his;fobj];
ts = [ts;etime(clock,initt)];

last_k = 100;
stop_accuracy = accuracy;
changes = ones(last_k,1)*inf;
fobj_old = fobj;


% L = norm(G)^2;
L = SpectralNorm(G);

for iter = 1:100000
    up = 0.5* norm(G*x-y)^2 + lambda * norm(x,1);
    [down,subgrad] =  topksum(x,k);
    down = down * lambda;
    subgrad = subgrad * lambda;
    fobj = up / down;
    
    % Solve the following optimization problem:
    % min_x 0.5 ||Gx-y||_2^2  + lambda ||x||_1 - fobj*<x-xt,subgrad>
    for in = 1:10
        grad = G'*(G*x-y) - fobj*subgrad;
        x = prox_l1(x-grad/L,lambda/L);
        if etime(clock,initt) > MaxTime,
            break;
        end
    end
    
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
