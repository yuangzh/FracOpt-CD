function best_t = nonconvex_prox_frac_topk(a,b,c,x,i,k)
% min_t (0.5 a t^2 + b t + c) / || x+tei ||_{top-k}

% case 1: x(i)+ t ²»ÊôÓÚtop-k set
% min_t (0.5 a t^2 + b t + c) / || x+tei ||_{top-k}
% min_t (0.5 a t^2 + b t + c)
t1 = -b/a;

% case 1: x(i)+ t ÊôÓÚtop-k set
% min_t (0.5 a t^2 + b t + c) / || x+tei ||_{top-k}
xJ = x;
xJ(i) = [];
delta = topksum(xJ,k-1);
% min_t (0.5 a t^2 + b t + c) / (delta + |x_i+t|) 

% case 1:A:
% x_i + t = 0
t2 = -x(i);
% case 1:B:
% x_i + t > 0
% min_t (0.5 a t^2 + b t + c) / (delta + x_i+t)
e = delta + x(i);
t3 = GetCriticalPoint(a,b,c,1,e);

% case 1:C:
% x_i + t < 0
% min_t (0.5 a t^2 + b t + c) / (delta - x_i-t) 
e = delta - x(i);
t4 = GetCriticalPoint(a,b,c,-1,e);


n = length(x);
absx = abs(x);
[sort_absx,ind]=sort(absx,'descend');



ts = [t1;t2;t3;t4;x;-x;0;x-x(i);-x-x(i)];




HandleObj = @(t)ComputeObj(t,a,b,c,x,i,k);
best_fobj = 1e100;
for i = 1:length(ts)
    test_obj = HandleObj(ts(i));
    if(test_obj<best_fobj)
        best_fobj = test_obj;
        best_t = ts(i);
    end
end






function [x] = GetCriticalPoint(a,b,c,d,e)
% (0.5 a tt + bt + c) / (dt + e)
% (at + b)(dt+e) = (0.5 a tt + bt + c) d
% ad tt + ate + bdt + be = 0.5 ad tt + bdt + cd
% ad tt + ate + be = 0.5 ad tt + cd
% 0.5 ad tt + aet  + be - cd = 0

A = 0.5*a*d;
B = a*e;
C = b*e - c*d;
t1 = (-B + sqrt(B*B - 4*A*C) )/(2*A);
t2 = (-B - sqrt(B*B - 4*A*C) )/(2*A);

A  =0.5*a;
B  = b;
C = c;
t3 = (-B + sqrt(B*B - 4*A*C) )/(2*A);
t4 = (-B - sqrt(B*B - 4*A*C) )/(2*A);

t1 = real(t1);
t2 = real(t2);
t3 = real(t3);
t4 = real(t4);
x = [t1;t2;t3;t4;-e/d];


function [f] = ComputeObj(t,a,b,c,x,i,k)
ei = 0*x;
ei(i) = 1;
r = x+t*ei;
[val,ind]=sort(abs(r),'descend');
% B = maxk(abs(r),k);
% sum(val(1:k)) - B
f = (0.5*a*t^2 + b*t + c) /  sum(val(1:k));

