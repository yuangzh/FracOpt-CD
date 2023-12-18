

function best_t = nonconvex_prox_frac_l2(a,b,c,x,i)
% min_t (0.5 a t^2 + b t + c) / || x+tei ||_2

% grad_down = (x_i + t) / norm(x+tei) 

% (at + b) o || x+tei ||_2 = (0.5 a t^2 + b t + c) o (x_i + t) / norm(x+tei)

% (at + b) o || x+tei ||^2_2 = (0.5 a t^2 + b t + c) o (x_i + t)

% (at + b) o (x'x + 2tx(i) + tt) = (0.5 a t^2 + b t + c) o (x_i + t)

% at o (x'x + 2tx(i) + tt) + b o (x'x + 2tx(i) + tt) = (0.5 a t^2 + b t + c) o x_i + (0.5 a t^2 + b t + c) o t

% at o (x'x + 2tx(i) + tt) + b o (x'x + 2tx(i) + tt) = (0.5 a t^2 + b t + c) o x_i + (0.5 a t^2 + b t + c) o t

xx = x'*x;
cof_3 = 0.5*a ;
cof_2 = 1.5*a*x(i);
cof_1 = a*xx + b*x(i) - c;
cof_0 = b*xx - c*x(i);
ts = roots([cof_3,cof_2,cof_1,cof_0]);
ts = real(ts);
ts = [ts;-x(i)];

HandleObj = @(t)(0.5*a*t^2 + b*t + c) /  sqrt(xx + 2*x(i)*t + t*t);
best_fobj = 1e100;
for i = 1:length(ts)
    test_obj = HandleObj(ts(i));
    if(test_obj<best_fobj)
        best_fobj = test_obj;
        best_t = ts(i);
    end
end







 
