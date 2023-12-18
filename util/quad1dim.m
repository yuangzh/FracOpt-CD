function s = quad1dim(u1,u2,lb,ub)
% min_{s}  u1 s + u2 s^2, s.t. lb <= s <= ub
% grad = u1 + 2 u2 * s = 0
% s = - u1 / (2*u2)

s = - u1 / (2*u2);
s = max(lb,min(s,ub));




