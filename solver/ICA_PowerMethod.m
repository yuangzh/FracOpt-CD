function [x,his,ts] = ICA_PowerMethod(A,x,speak,accuracy,MaxTime)
% max_x ||Ax||_4^4, s.t. ||x||_2 = 1
initt = clock;
x = x / norm(x);
his = [];




his = []; ts =[];

HandleObj = @(x)(x'*x) / (norm(A*x,4)^2);
fobj = HandleObj(x);
his = [his;fobj];
ts = [ts;etime(clock,initt)];

last_k = 100;
stop_accuracy = accuracy;
changes = ones(last_k,1)*inf;
fobj_old = fobj;



for iter = 1:500000
    [down_4,grad] = ComputeObjGradNorm4(x,A);
    fobj = (x'*x) / sqrt(down_4);
    
    % min_x <x,-grad>, s.t. ||x|| = 2
    % min_x <x,-grad> + 0.5 x'x, s.t. ||x|| = 2
    % min_x 0.5 ||x - grad||_2^2, s.t. ||x|| = 2
 
    x = grad;
    x = x / norm(x);
%     if(norm(x-xt)<acc),break;end
    if(speak),fprintf('iter:%d, fobj:%f\n',iter,fobj);end
    
    
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

