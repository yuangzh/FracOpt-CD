function best_t = nonconvex_prox_frac_l1(a,b,c,g,d)

HandleObj = @(t) (0.5*a*t*t + b*t + c) / (norm(g*t+d,1)^2);



I = (g==0);
g(I)=[];d(I)=[];
m = length(g);

v = abs(d)./abs(g);
v = sort(v);
v = [-1e10;-flip(v);v;1e10];
lenv = length(v);


his= [];
ts = [];


 
for i = 1:(lenv-1)
    tt = (v(i)+v(i+1)) / 2;
    o  = sign(g*tt + d);
    gamma = o'*d;
    phi = o'*g;
    % (a t + b) * (phi t + gamma) = (a t t + 2bt + 2c) phi
    % phi t a t + gamma a t + phi b t + gamma b = a phi t t + 2b phi t + 2cphi 
    % gamma a t + phi b t - 2b phi t =  2cphi -  gamma b
    % t =  (2cphi -  gamma b) / (gamma a  + phi b  - 2b phi)
    t = (2 * c * phi - gamma * b) / (gamma*a  + phi*b  - 2*b*phi);
   if(t>v(i) && t<v(i+1))
    ts = [ts;t];
    his = [his;HandleObj(t)];
   end
end


ss = [d./g;-d./g];
for i=1:length(ss)
    ts = [ts;ss(i)];
    his = [his;HandleObj(ss(i))];
end


[~,ind]=min(his);
best_t = ts(ind);


