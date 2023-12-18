function L = SpectralNorm(G)

[m,n] = size(G);
if(n<m)
    [~,L] = laneig(G'*G,1,'AL');
else
    [~,L] = laneig(G*G',1,'AL');
end