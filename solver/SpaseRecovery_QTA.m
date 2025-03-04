function [x,his,ts] = SpaseRecovery_QTA(G,y,x,lambda,k,speak,accuracy,MaxTime)
% min_x (0.5 ||Gx-y||_2^2  + lambda ||x||_1)  /  (lambda ||x||_{top_k})
% min_{x,t} t^2 (0.5 ||Gx-y||_2^2  + lambda ||x||_1) - 2t sqrt (lambda ||x||_{top_k})

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
for iter = 1:10000000
    up = 0.5* norm(G*x-y)^2 + lambda * norm(x,1);
    [down1,subgrad1] =  topksum(x,k);
    down = down1 * lambda;
    subgrad = subgrad1 * lambda;
    fobj = up / down;
    
    beta = sqrt(down) / up;

    % Solve the following optimization problem:
    % min_x beta^2 (0.5 ||Gx-y||_2^2  + lambda ||x||_1) - 2*beta sqrt(lambda ||x||_{top_k})
    % min_x  (0.5 ||Gx-y||_2^2  + lambda ||x||_1) - 2/beta sqrt(lambda ||x||_{top_k})
    zzz = -2/beta*0.5* (down)^(-0.5)* subgrad;

    grad = G'*(G*x-y) + zzz;
    x = prox_l1(x-grad/L,lambda/L);

   
 
%     HandleObjSmooth = @(x)computeObj(x,G,y,zzz);
%     HandleObjNonSmooth = @(x)lambda*sum(abs(x));
%     HandleProx = @(theta,a)computeprox(theta,a,lambda);
%     max_iter = 1;
%     Lip = 1;
%     x = ProximalGradient(x,HandleObjSmooth,HandleObjNonSmooth,HandleProx,max_iter,Lip);
    
    
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




function [fobj,grad] = computeObj(x,G,y,zzz)
diff = G*x-y;
fobj = 1/2*norm(diff)^2 + x'*zzz;
grad = G'*diff + zzz;

function [x] = computeprox(theta,a,lambda)
% 0.5 theta ||x - a||^2 + g(x)
[x] = threadholding_l1(a,lambda/theta);
% x = prox_l1(a-grad/theta,lambda/theta);

function [x] = threadholding_l1(a,lambda)
% solving the following OP:
% min_{x} 0.5 ||x - a||^2 + lambda * sum(abs(x))
x = sign(a).*max(0,abs(a)-lambda);
