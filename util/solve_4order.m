function [t_best] = solve_4order(cof1,cof2,cof3,cof4)
% min_t cof_1*t^1 + cof_2*t^2 + cof_3*t^3 + cof_4*t^4


cof5 = 1e-5;
HandleObj = @(t)cof1.*t + cof2.*t.^2 + cof3.*t.^3 + cof4.*t.^4 + cof5.*abs(t).^5;
% t1_best = fminbnd(HandleObj,-R,R);

% grad_t = cof_1 + 2*cof_2*t + 3*cof_3*t^2 + 4*cof_4*t^3 = 0
ts1 = roots([5*cof5,4*cof4,3*cof3,2*cof2,cof1]);
ts2 = roots([-5*cof5,4*cof4,3*cof3,2*cof2,cof1]);
ts = [real(ts1);real(ts2)];
fs = HandleObj(ts);
[~,ind]=min(fs);
t_best = ts(ind);