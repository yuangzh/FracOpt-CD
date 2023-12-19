function x_best = nonconvex_prox_frac_l2l42222(a2,a1,a0,b4,b3,b2,b1,b0,L,U)
%                            (a2 t^2 + a1 t + a0) + C t^2
% min_{L<=t<=U} ---------------------------------------------
%               £¨b4 t^4 + b3 t^3 + b2 t^2 + b1 t^1 + b0£©
% We assume that the denominator is positive
% a2,a1,a0,b4,b3,b2,b1,b0,L,U

Proj = @(x) min(U,max(L,x));
HandleObj = @(t)ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0);
cof_5 = 4*a2.*b4;
cof_4 = 4*a2.*b3 + 2*a1.* b4;
cof_3 = 4*a2.*b2 + 2*a1.*b3; 
cof_2 = 4*a2.*b1+ (2*a1) .* ( b2);
cof_1 = 4*a2.*b0  + (2*a1) .* (b1);
cof_0 = 2*a1.*b0;
% left = cof_5 * t^5 + cof_4*t^4 + cof_3*t^3 + cof_2*t^2 +  cof_1*t   + cof_0;




cof2_5 = (a2) .* ( 4* b4 );
cof2_4 = (a2) .* (3* b3) +    ( a1 ) .* ( 4* b4);
cof2_3 = (a2) .* (2 *b2 ) +    ( a1) .* ( 3* b3 )  +   ( a0) .* ( 4* b4 ); 
cof2_2 = (a2 ) .*  b1 +    ( a1) .* ( 2 *b2 )  +   ( a0) .* ( 3* b3);
cof2_1 =   ( a1) .* ( b1)  +   ( a0) .* ( 2 *b2 );
cof2_0 =   ( a0) .* ( b1);

Left = @(t)(4.*a2.*t + 2.*a1) .* (b4 *t.^4 + b3.* t.^3 + b2.*t.^2 + b1.*t + b0);
Right = @(t)(a2.*t.^2 + a1.*t + a0) .* ( 4.* b4.* t.^3 + 3.*b3.*t.^2 + 2 *b2.*t + b1);

c4 = cof_4-cof2_4; c3 = cof_3-cof2_3; c2 = cof_2-cof2_2; c1 = cof_1-cof2_1; c0 = cof_0-cof2_0;


% find t such that left == right
%  cof_5 * t^5 + cof_4*t^4 + cof_3*t^3 + cof_2*t^2 +  cof_1*t   + cof_0  = cof2_5 * t^5 + cof2_4*t^4 + cof2_3*t^3 + cof2_2*t^2 +  cof2_1*t   + cof2_0
% it reduces to 
%  cof_4*t^4 + cof_3*t^3 + cof_2*t^2 +  cof_1*t   + cof_0  = cof2_4*t^4 + cof2_3*t^3 + cof2_2*t^2 +  cof2_1*t   + cof2_0

xs1 = roots([c4;c3;c2;c1;c0]);
xs2 = roots([b4;b3;b2;b1;b0]);

xs1 = [abs(xs1);real(xs1)];
xs2 = [abs(xs2);real(xs2)];

xs = [0;xs1;xs2];


xs = [L;U; Proj(xs)];

for iii=1:length(xs)
fs(iii) = HandleObj(xs(iii));
end

[~,ind] = min(fs);
x_best = xs(ind);

% ind = find(fs==fs(ind));
% x_candidate = xs(ind);
% [~,ind2] = min(abs(x_candidate));
% x_best = x_candidate(ind2);





function [fobj,grad] = ComputeObj(t,a2,a1,a0,b4,b3,b2,b1,b0)
up =  a2*t*t + a1*t + a0;
ooo = b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0;
if(ooo<0)
    fobj = 1e100;
    grad = 1e100;
    return;
end
down = sqrt(ooo);
fobj = up./down;
grad_up = 2*a2*t + a1;
grad_down = 0.5*(b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0)^(-0.5)*(4*b4*t*t*t + 3*b3*t*t + 2*b2*t + b1);
grad = grad_up - fobj.*grad_down;
grad = grad ./ down;

% 2*a2*t + a1 = (a2*t*t + a1*t + a0) / sqrt(b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0)  o  0.5*(b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0)^(-0.5)*(4*b4*t*t*t + 3*b3*t*t + 2*b2*t + b1)
% 2*(2*a2*t + a1) (b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0) = (a2*t*t + a1*t + a0) (4*b4*t*t*t + 3*b3*t*t + 2*b2*t + b1)
% 2*(2*a2*t + a1) (b4*t*t*t*t + b3*t*t*t + b2*t*t + b1*t + b0) = (a2*t*t + a1*t + a0) (4*b4*t*t*t + 3*b3*t*t + 2*b2*t + b1)
