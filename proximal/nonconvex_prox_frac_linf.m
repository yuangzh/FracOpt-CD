function best = nonconvex_prox_frac_linf(a,b,c,g,d)
% min_t  (a t t + b t + c ) / (norm(g t + d,inf)^2)
% We assume that (at^2 + bt + c) >=0 for all t, otherwise 
% this problem cannot be solved exactly usning breakpoint search method


p = [g;-g]; w = [d;-d]; 
clear g d;
% min_t  (a t t + b t + c ) / (max(pt+w)^2)
HandleObj = @(t) (a*t*t + b*t + c) / (max(p*t+w)^2);
% best2 = fminsearch(HandleObj,0);

his= [];
ts = [];
for i = 1:length(w),
    ttt   =  - (b * w(i)  - c*2*p(i))  /   ( 2*a*w(i) - b*p(i));
    ts = [ts;ttt];
    his = [his;HandleObj(ttt)];
end

 
[~,ind]=min(his);
best = ts(ind);
% fff2 = HandleObj(best2);
% fff1 = HandleObj(best);
% if(fff2+1e-4<fff1)
%     fff1
%     fff2
%  
%     ddd
% end
function [fobj,grad] = ComputeObj_t(t,a,b,c,p,w)
% (a t t + b t + c ) / (max(pt+w)^2)
o = p*t+w;
[val,ind] = max(o);
up = a*t*t + b*t + c;
down = (val^2);
grad_up = 2*a*t + b;
grad_down = 2*val*p(ind);
fobj = up / down;
grad = (grad_up * down - up*grad_down) / (down^2);

grad2 = grad_up / down - up*grad_down / (down^2);

 
grad =  ( 2*a*w(ind) - b*2*p(ind) + b * p(ind) ) *t +        b * w(ind)  - c*2*p(ind)


% (2*a*t + b) / (val^2) = (a*t*t + b*t + c) * 2*val*p(ind) / (val^4)
% (2*a*t + b) * val = (a*t*t + b*t + c) * 2*p(ind) 
% (2*a*t + b) * (p_i t+w_i) = (a*t*t + b*t + c) * 2*p(i) 





 
function [s] = quadfrac2(u0,u1,u2,d0,d1,d2)
% min_{s}  (u0 + u1 s + u2 s^2)  /  (d0 + d1s + d2 s^2)
Handle = @(s) (u2.*s.*s + u1.*s + u0)  ./  (d2.*s.*s + d1.*s + d0);
% grad_s =  (2 u2 s+u1)*(d2 s^2 + d1 s + d0) - (u2 s^2 + u1 s + u0)*(2 d2 s + d1) =0
% (2 u2 s + u1)*(d2 s^2 + d1 s + d0) - (2 d2 s + d1)*(u2 s^2 + u1 s + u0) =0
% s = fminsearch(Handle,0);
% s
% Handle(s)
% return;
c2 = u2.*d1 - d2.*u1;
c1 = 2.*u2.*d0 - 2.*d2.*u0;
c0 = u1.*d0 - d1.*u0;

% c2 x^2 + c1 x + c0 = 0
dd = sqrt(c1.*c1 - 4*c2.*c0);
dd = real(dd);
x1 = (-c1+dd)./(2*c2);
x2 = (-c1-dd)./(2*c2);
ind = find(c2==0);
x1(ind) = -c0(ind)./c1(ind);
x2(ind) = x1(ind);
nanindex = find(isnan(x1) | isinf(x1));
x1(nanindex)=0; x2(nanindex)=0;
xx = [x1 x2]';
fs = [Handle(x1) Handle(x2)];
[fs,ind] = min(fs');
s = xx([0:2:(2*length(u0)-1)]+ind)';
fs = fs(:);



