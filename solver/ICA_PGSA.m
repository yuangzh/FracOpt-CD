function [x,his,ts] = ICA_PGSA(A,x,speak,accuracy,MaxTime)
% min_x x'x / sqrt(||Ax||_4^4)

initt = clock;
his = []; ts =[];
fobj = (x'*x) / sqrt(ComputeObjGradNorm4(x,A));
his = [his;fobj];
ts = [ts;etime(clock,initt)];

last_k = 100;
stop_accuracy = accuracy;
changes = ones(last_k,1)*inf;
fobj_old = fobj;


for iter = 1:50000
    [fobj1,grad] = ComputeObjGradNorm4_sqrt(x,A);
    fobj = (x'*x) / fobj1;
    % solve the following problem:
    % min_x x'x - 2*0.5*lambda g'(x-xt)
    % min_x x'x - 2*(0.5*lambda g)'x
    x = 0.5*fobj*grad;
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
