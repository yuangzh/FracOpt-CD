function [x,his,ts] = ICA_FractionalCD(A,x,speak,accuracy,MaxTime)
% min_x x'x / sqrt(||Ax||_4^4)
% min_t (x+tei)'(x+tei) / sqrt(||A(x+tei)||_4^4)

initt = clock;
HandleObj = @(x)(x'*x) / (norm(A*x,4)^2);

theta = 1e-5;




[m,n] = size(A);


his = []; ts =[];
fobj = HandleObj(x);
his = [his;fobj];
ts = [ts;etime(clock,initt)];

last_k = 5000;
stop_accuracy = accuracy;
changes = ones(last_k,1)*inf;
fobj_old = fobj;

Ax = A*x;

for iter = 1:1000000

    fobj =  (x'*x) / (norm(Ax,4)^2);
    his(iter) = fobj;

    % solve the following problem 
    % i = randperm(n,1);
    i = mod(iter,n)+1;

a0 = x'*x;
a1 = 2*x(i);
a2 = 1;

 
  [b0,b1,b2,b3,b4] = ComputeCofPoly4(x,A,i,Ax);
 

  alpha = nonconvex_prox_frac_l2l4(a2+theta,a1,a0,b4,b3,b2,b1,b0,-1e10,1e10);
%      ei = zeros(n,1);
%     ei(i) = 1;
%    HandleObj_alpha  = @(alpha) (norm(x+alpha*ei)^2) / ( norm(A*(x+alpha*ei),4)^2 ) ;
%       [alpha,~,~,~] = fminsearch(HandleObj_alpha,0); 

% (a2 t^2 + a1 t + a0)^2
% --------------------------------------
% b4 t^4 + b3 t^3 + b2 t^2 + b1 t^1 + b0

%       if(HandleObj_alpha(alpha) > HandleObj_alpha(0) + 0.01)
%           HandleObj_alpha(alpha)
%           HandleObj_alpha(0)
%           ddd
%       end

  x(i) = x(i) + alpha;
  % update Ax
  Ax = Ax + alpha*A(:,i);
 
   if(speak),
         fprintf('iter:%d, fobj:%f, alpha:%f\n',iter,fobj,alpha);
   end

    his(iter) = fobj;
    ts(iter) = etime(clock,initt);
    
    if etime(clock,initt) > MaxTime,
        break;
    end
    if(fobj>fobj_old+1e-3)
        dddd
    end

    rel_change = abs((fobj - fobj_old)/max(1,fobj_old));
    changes = [changes(2:end);rel_change];
    fobj_old = fobj;
    if(mean(changes)<stop_accuracy),break;end
     
end
